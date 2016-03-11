//
// Created by ryan on 3/7/16.
//

#ifndef EXPLICIT_FEA_SETUP_H
#define EXPLICIT_FEA_SETUP_H

#include <memory>
#include <rapidjson/document.h>
#include <rapidjson/filereadstream.h>
#include "explicit_system.h"

namespace explicit_fea {
    /**
     * Opens the specified json file and parses the data into a `rapidjson::Document` and returns the result.
     * The config document should have key's "nodes", "elems", "props" and "bcs". Optionally, there can be keys
     * "forces" for prescribed forces, "t_0" for the initial time of the system,
     * "initial_displacements" and "initial_velocities" for initial conditions (assumed to be zero if not provided),
     * and an "options" key specifying additional paramters to the system.
     *
     * @param config_filename `std::string`. The location of the configuration json file.
     * @return Document`rapidjson::Document`
     */
    rapidjson::Document parseJSONConfig(const std::string &config_filename);

    /**
     * Parses the file indicated by the "nodes" key in `config_doc` into a vector of `explicit_fea::Node`'s.
     *
     * @param config_doc `rapidjson::Document`. Document storing the file name containing the nodal coordinates.
     * @return Nodal coordinates. `std::vector<explicit_fea::Node>`.
     */
    std::vector<Node> createNodeVecFromJSON(const rapidjson::Document &config_doc);

    /**
     * Parses the files indicated by the "elems" and "props" keys in `config_doc` into a
     * vector of `explicit_fea::Elem` pointers.
     *
     * @param config_doc `rapidjson::Document`. Document storing the file names of the csv files that contain
     *                    the node number designations for each element and elemental properties.
     * @return Elements. `std::vector<explicit_fea::BeamElement>`.
     */
    std::vector<std::unique_ptr<BeamElement>> createElemVecFromJSON(const rapidjson::Document &config_doc);

    /**
     * Parses the file indicated by the "bcs" key in `config_doc` into a vector of `explicit_fea::BC` pointers.
     *
     * @param config_doc `rapidjson::Document`. Document storing the file name containing the boundary conditions.
     * @return Boundary conditions. `BCList`.
     */
    BCList createBCVecFromJSON(const rapidjson::Document &config_doc);

    /**
     * Parses the file indicated by the "forces" key in `config_doc` into a vector of `explicit_fea::Forces` pointers.
     *
     * @param config_doc `rapidjson::Document`. Document storing the file name containing the prescribed forces.
     * @return External forces. `ForceList`.
     */
    ForceList createForceVecFromJSON(const rapidjson::Document &config_doc);

    /**
     * @brief Parses the file indicated by the specified key into a `ColumnVector`.
     * @details The parsed data must have 1 value per line and the number of lines must equal the specified `size`.
     *
     * @param config_doc Document storing a filename associated with the given `key`.
     * @param key Key in the document that holds the filename to load.
     * @param size The expected size of the parsed vector.
     * @return Column vector, used internally to load initial conditions from CSV files.
     */
    ColumnVector createColumnVectorFromJSON(const rapidjson::Document &config_doc, std::string key, unsigned int size);
}

#endif //EXPLICIT_FEA_SETUP_H
