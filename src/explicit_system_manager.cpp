#include "explicit_system_manager.h"
#include <boost/format.hpp>
#include <chrono>
#include <rapidjson/ostreamwrapper.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "setup.h"

namespace explicit_fea {
    ExplicitSystemManager::ExplicitSystemManager(const std::string& _config_filename)
    {
        config_doc = parseJSONConfig(_config_filename);
        options.load(config_doc);
        if (config_doc.HasMember("iteration_number"))
            iteration_number = config_doc["iteration_number"].GetUint();
        else {
            iteration_number = 0;
            rapidjson::Value i;
            i.SetUint(iteration_number);
            config_doc.AddMember("iteration_number", i, config_doc.GetAllocator());
        }
        init_time_from_config();
        construct_system();
    }

    void ExplicitSystemManager::Options::load(const rapidjson::Document &config_doc) {
        if (config_doc.HasMember("options")) {
            if (config_doc["options"].HasMember("state_filename")) {
                if (!config_doc["options"]["state_filename"].IsString()) {
                    throw std::runtime_error("state_filename provided in options configuration is not a string.");
                }
                state_filename = config_doc["options"]["state_filename"].GetString();
            }
            if (config_doc["options"].HasMember("nodal_displacements_filename")) {
                if (!config_doc["options"]["nodal_displacements_filename"].IsString()) {
                    throw std::runtime_error("nodal_displacements_filename provided in options configuration is not a string.");
                }
                nodal_displacements_filename = config_doc["options"]["nodal_displacements_filename"].GetString();
            }
            if (config_doc["options"].HasMember("nodal_velocities_filename")) {
                if (!config_doc["options"]["nodal_velocities_filename"].IsString()) {
                    throw std::runtime_error("nodal_velocities_filename provided in options configuration is not a string.");
                }
                nodal_velocities_filename = config_doc["options"]["nodal_velocities_filename"].GetString();
            }
            if (config_doc["options"].HasMember("save_frequency")) {
                if (!config_doc["options"]["save_frequency"].IsNumber()) {
                    throw std::runtime_error("save_frequency provided in options configuration is not a number.");
                }
                save_frequency = config_doc["options"]["save_frequency"].GetUint();
            }
            if (config_doc["options"].HasMember("verbose")) {
                if (!config_doc["options"]["verbose"].IsBool()) {
                    throw std::runtime_error("verbose provided in options configuration is not a bool.");
                }
                verbose = config_doc["options"]["verbose"].GetBool();
            }
        }
    }

	ExplicitSystem const& ExplicitSystemManager::getExplicitSystem() const {
		return *explicit_system;
	}

	rapidjson::Document const& ExplicitSystemManager::getConfigDoc() const {
		return config_doc;
	}

	double ExplicitSystemManager::getTimeStep() const {
		return dt;
	}

	unsigned int ExplicitSystemManager::getIterationNumber() const {
		return iteration_number;
	}

    void ExplicitSystemManager::run() {
        if (options.verbose)
            std::cout << "Starting analysis:" << std::endl;
		auto t1 = std::chrono::high_resolution_clock::now();
        ValueCompare<double> compare;
        double time_period = end_time - start_time;
        int old_percent = 0;
        int new_percent = 0;
        if (options.verbose) {
            std::cout << "\tSaving initial system state..." << std::endl;
        }
        dump_system();

        if (options.verbose)
            std::cout << "\tAdvancing equations of motion..." << std::endl;
        while (compare.lessThan(explicit_system->getTime(), end_time)) {
            explicit_system->update(dt);
            ++iteration_number;

            if (options.save_frequency > 0 && iteration_number % options.save_frequency == 0)
                dump_system();

			if (options.verbose) {
                new_percent = (int)((explicit_system->getTime() - start_time) / time_period * 100 + 0.1);
                if (new_percent > old_percent) {
                    std::cout << "\t\t" << new_percent << "% completed." << std::endl;
                    old_percent = new_percent;
                }
			}
        }
        if (options.verbose) {
            std::cout << "\tSaving final system state..." << std::endl;
        }
        dump_system();
		auto t2 = std::chrono::high_resolution_clock::now();
		auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
		if (options.verbose)
			std::cout << "Explicit time integration completed in " << time_span.count() << " seconds." << std::endl;
    }

    void ExplicitSystemManager::init_time_from_config() {
        if (config_doc.HasMember("start_time")) {
            if (!config_doc["start_time"].IsNumber()) {
                throw std::runtime_error("start_time provided in options configuration is not a number.");
            }
            start_time = config_doc["start_time"].GetDouble();
        }
        else {
            throw std::runtime_error("Configuration file does not have requested member variable start_time.");
        }
        if (config_doc.HasMember("end_time")) {
            if (!config_doc["end_time"].IsNumber()) {
                throw std::runtime_error("end_time provided in options configuration is not a number.");
            }
            end_time = config_doc["end_time"].GetDouble();
        }
        else {
            throw std::runtime_error("Configuration file does not have requested member variable end_time.");
        }
    }

    void ExplicitSystemManager::dump_system() {
        int denominator = std::max((int)options.save_frequency, 1);
        std::string filename_tail = (boost::format("_%05d") % (iteration_number / denominator)).str();
        std::string displacements_filename = (boost::format("%s%s.txt") % options.nodal_displacements_filename % filename_tail).str();
        std::string velocities_filename = (boost::format("%s%s.txt") % options.nodal_velocities_filename % filename_tail).str();
		std::string forces_filename = (boost::format("%s%s.txt") % options.nodal_forces_filename % filename_tail).str();

        save_column_vector(explicit_system->getDisplacements(), displacements_filename);
        config_doc["nodal_displacements"].SetString(displacements_filename.c_str(), displacements_filename.length());

        save_column_vector(explicit_system->getVelocities(), velocities_filename);
        config_doc["nodal_velocities"].SetString(velocities_filename.c_str(), velocities_filename.length());

		save_column_vector(explicit_system->getForces(), forces_filename);

        config_doc["start_time"].SetDouble(explicit_system->getTime());
        config_doc["iteration_number"].SetUint(iteration_number);

        std::string state_filename((boost::format("%s%s.json") % options.state_filename % filename_tail).str());
        std::ofstream out_file(state_filename);
        if (out_file.is_open()) {
            rapidjson::OStreamWrapper o_wrapper(out_file);
            rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(o_wrapper);
            config_doc.Accept(writer);
            out_file.close();
        }
        else {
            throw std::runtime_error((boost::format("Unable to open %s.") % state_filename).str());
        }
    }

    void ExplicitSystemManager::save_column_vector(const ColumnVector &vec, std::string filename) {
        std::ofstream out_file(filename);
        if (out_file.is_open()) {
            for (int i = 0; i < vec.size(); ++i) {
                out_file << std::setprecision(15) << vec[i] << "\n";
            }
            out_file.close();
        }
        else {
            throw std::runtime_error((boost::format("Unable to open %s.") % filename).str());
        }
    }

    void ExplicitSystemManager::construct_system() {
        if (options.verbose) {
            std::cout << "Beginning construction of system:" << std::endl;
        }
		auto t1 = std::chrono::high_resolution_clock::now();
        if (options.verbose) {
            std::cout << "\tParsing node list..." << std::endl;
        }
        std::vector<Node> nodes = createNodeVecFromJSON(config_doc);
        if (options.verbose) {
            std::cout << "\tParsing element list..." << std::endl;
        }
        std::vector<std::unique_ptr<BeamElement>> elems = createElemVecFromJSON(config_doc);
        if (options.verbose) {
            std::cout << "\tParsing node boundary conditions..." << std::endl;
        }
        BCList bcs = createBCVecFromJSON(config_doc);
        if (options.verbose) {
            std::cout << "\tCreating mesh..." << std::endl;
        }
        Mesh mesh(nodes, elems, std::move(bcs));
        if (options.verbose) {
            std::cout << "\tParsing external forces..." << std::endl;
        }
        ForceList forces = createForceVecFromJSON(config_doc);
        if (options.verbose) {
            std::cout << "\tParsing initial conditions..." << std::endl;
        }
        ColumnVector initial_displacements = createColumnVectorFromJSON(config_doc,
                                                                        "nodal_displacements",
                                                                        mesh.get_global_stiffness_matrix().cols());
        ColumnVector initial_velocities = createColumnVectorFromJSON(config_doc,
                                                                     "nodal_velocities",
                                                                     mesh.get_global_stiffness_matrix().cols());
        if (options.verbose) {
            std::cout << "\tMoving data into system..." << std::endl;
        }
        ExplicitSystem::Options es_options;
        options.load(config_doc);
        std::unique_ptr<ExplicitSystem> temp_es =
                std::unique_ptr<ExplicitSystem>(new ExplicitSystem(std::move(mesh), std::move(forces),
                                                                   initial_displacements, initial_velocities,
                                                                   start_time, es_options));
        std::swap(temp_es, explicit_system);

        dt = estimate_stable_timestep(nodes, elems);
		auto t2 = std::chrono::high_resolution_clock::now();
		auto time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

		if (options.verbose) {
			std::cout << "Constructed system of " << elems.size() <<" elements and " << nodes.size() << " nodes in " << time_span.count() << " seconds.\n" 
			          << "Estimated stable time step is " << dt << " seconds." << std::endl;
		}
    }

} // namespace explicit_fea

