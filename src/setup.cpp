//
// Created by ryan on 3/7/16.
//

#include "setup.h"
#include <boost/format.hpp>
#include <rapidjson/error/en.h>
#include "csv_parser.h"
#include "timoshenko_beam_element.h"

namespace explicit_fea {
    namespace {
        template<typename T>
        void createVectorFromJSON(const rapidjson::Document &config_doc,
                                  const std::string &variable,
                                  std::vector<T> &data) {
            if (!config_doc.HasMember(variable.c_str())) {
                throw std::runtime_error(
                        (boost::format("Configuration file does not have requested member variable %s.") %
                         variable).str()
                );
            }
            if (!config_doc[variable.c_str()].IsString()) {
                throw std::runtime_error(
                        (boost::format("Value associated with variable %s is not a string.") % variable).str()
                );
            }
            CSVParser csv;
            std::string filename(config_doc[variable.c_str()].GetString());
            csv.parseToVector(filename, data);
            if (data.size() == 0) {
                throw std::runtime_error(
                        (boost::format("No data was loaded for variable %s.") % variable).str()
                );
            }
        }
    }

    rapidjson::Document parseJSONConfig(const std::string &config_filename) {
        rapidjson::Document config_doc;

        FILE *config_file_ptr = fopen(config_filename.c_str(), "r");

        if (!config_file_ptr) {
            throw std::runtime_error(
                    (boost::format("Cannot open configuration input file %s.") % config_filename).str()
            );
        }
        char readBuffer[65536];
        rapidjson::FileReadStream config_stream(config_file_ptr, readBuffer, sizeof(readBuffer));
        config_doc.ParseStream(config_stream);
        if (config_doc.HasParseError()) {
            throw std::runtime_error(
                    (boost::format("Error parsing %s (offset %d):\t%s") % config_filename % config_doc.GetErrorOffset() % rapidjson::GetParseError_En(config_doc.GetParseError())).str()
            );
        }
        fclose(config_file_ptr);
        return config_doc;
    }

    std::vector<Node> createNodeVecFromJSON(const rapidjson::Document &config_doc) {
        std::vector<std::vector<double> > nodes_vec;
        createVectorFromJSON(config_doc, "nodes", nodes_vec);
        std::vector<Node> nodes_out(nodes_vec.size());
        Node n;

        for (size_t i = 0; i < nodes_vec.size(); ++i) {

            if (nodes_vec[i].size() != 3) {
                throw std::runtime_error(
                        (boost::format("Row %d in nodes does not specify x, y and z coordinates.") % i).str()
                );
            }
            n << nodes_vec[i][0], nodes_vec[i][1], nodes_vec[i][2];
            nodes_out[i] = n;
        }
        return nodes_out;
    }

    std::vector<std::unique_ptr<BeamElement>> createElemVecFromJSON(const rapidjson::Document &config_doc) {
        std::vector< std::vector<int> > elems_vec;
        std::vector< std::vector<double> > props_vec;

        createVectorFromJSON(config_doc, "elems", elems_vec);
        createVectorFromJSON(config_doc, "props", props_vec);

        if (elems_vec.size() != props_vec.size()) {
            throw std::runtime_error("The number of rows in elems did not match props.");
        }

        std::vector<std::unique_ptr<BeamElement>> elems_out(elems_vec.size());
        Eigen::Vector3d normal_vec;
        for (size_t i = 0; i < elems_vec.size(); ++i) {
            if (elems_vec[i].size() != 2) {
                throw std::runtime_error(
                        (boost::format("Row %d in elems does not specify 2 nodal indices [nn1,nn2].") % i).str()
                );
            }
            if (props_vec[i].size() != 10) {
                throw std::runtime_error(
                        (boost::format("Row %d  in props does not specify the 10 property values "
                                       "[youngs_modulus,shear_modulus,area,Iz,Iy,J,density,nx,ny,nz]") % i).str()
                );
            }
            normal_vec << props_vec[i][7], props_vec[i][8], props_vec[i][9];
            Props p(props_vec[i][0], props_vec[i][1], props_vec[i][2],
                    props_vec[i][3], props_vec[i][4], props_vec[i][5],
                    props_vec[i][6], normal_vec);

            elems_out[i] = std::unique_ptr<TimoshenkoBeamElement>(new TimoshenkoBeamElement(elems_vec[i][0], elems_vec[i][1], p));
        }
        return elems_out;
    }

    std::vector<std::unique_ptr<BC>> createBCVecFromJSON(const rapidjson::Document &config_doc) {
        std::vector< std::vector<double> > bcs_vec;

        createVectorFromJSON(config_doc, "bcs", bcs_vec);
        std::vector<std::unique_ptr<BC>> bcs_out(bcs_vec.size());

        for (size_t i = 0; i < bcs_vec.size(); ++i) {
            if (bcs_vec[i].size() != 4) {
                throw std::runtime_error(
                        (boost::format("Row %d in bcs does not specify [node number,DOF,value,type].") % i).str()
                );
            }
            bcs_out[i] = std::unique_ptr<ConstantBC>(new ConstantBC((unsigned int) bcs_vec[i][0], (unsigned int) bcs_vec[i][1], bcs_vec[i][2], (BC::Type)(unsigned int) bcs_vec[i][3]));
        }
        return bcs_out;
    }

    std::vector<std::unique_ptr<Force>> createForceVecFromJSON(const rapidjson::Document &config_doc) {
        std::vector<std::unique_ptr<Force>> forces_out;
        if (config_doc.HasMember("forces")) {
            std::vector<std::vector<double> > forces_vec;
            
            createVectorFromJSON(config_doc, "forces", forces_vec);

            forces_out.resize(forces_vec.size());
            for (size_t i = 0; i < forces_vec.size(); ++i) {
                if (forces_vec[i].size() != 3) {
                    throw std::runtime_error(
                            (boost::format("Row %d in forces does not specify [node number,DOF,value].") % i).str()
                    );
                }

                forces_out[i] = std::unique_ptr<ConstantForce>(new ConstantForce((unsigned int) forces_vec[i][0], (unsigned int) forces_vec[i][1], forces_vec[i][2]));
            }
        }
        return forces_out;
    }

    ColumnVector createColumnVectorFromJSON(const rapidjson::Document &config_doc, std::string key, unsigned int size) {
        ColumnVector col_vec(size);

        if (config_doc.HasMember(key.c_str())) {
            std::vector< std::vector<double> > vec;
            
            createVectorFromJSON(config_doc, key, vec);

            if (vec.size() != size) {
                throw std::runtime_error(
                        (boost::format("Key specified by %s does not have the required %d values. %d entries were parsed.")
                         % key % size % vec.size()).str()
                );
            }
            for (int i = 0; i < size; ++i) {
                if (vec[i].size() != 1) {
                    throw std::runtime_error(
                            (boost::format("Row %d in the file specified by %s has %d values. There should only be one nodal value per line.")
                             % i % key % size % vec[i].size()).str()
                    );
                }
                col_vec[i] = vec[i][0];
            }
        }
        else {
            col_vec.setZero();
        }

        return col_vec;
    }
}
