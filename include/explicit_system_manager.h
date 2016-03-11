

#ifndef EXPLICIT_FEA_EXPLICIT_SYSTEM_MANAGER_H
#define EXPLICIT_FEA_EXPLICIT_SYSTEM_MANAGER_H

#include "explicit_system.h"

namespace explicit_fea {

    /**
     * @brief Manages explicit time integration
     * @details This class is responsible for parsing a configuration
     * file, constructing the system, and running the analysis.
     */
    class ExplicitSystemManager {
    public:
        struct Options {
            Options() :
                    verbose(false),
                    save_frequency(0),
                    state_filename("state"),
                    nodal_displacements_filename("nodal_displacements"),
                    nodal_velocities_filename("nodal_velocities"),
					nodal_forces_filename("nodal_forces")

            {}

            /**
             * @brief Updates member variables based on available keys in the DOM.
             * @details If keys are present in the DOM that match member variable names
             * the value is updated to the value contained in the document.
             *
             * @param config_doc DOM containing optional values.
             */
            void load(const rapidjson::Document &config_doc);

            /**
             * @brief Flag to print status updates during simulation, default = `false`.
             * @details If `True` status updates will be written to the console at the frequency
             * given by `save_frequency`.
             */
            bool verbose;

            /**
             * @brief Frequency at which to save the state of the model, default = `0`.
             * @details The state of the system will be saved every time
             * `i % save_frequency == 0` where `i` is the current iteration of the system.
             * Set to 0 to disable saving intermediate system states during simulation.
             *
             */
            unsigned int save_frequency;

            /**
             * @brief Prefix to filename when saving system state, default = `"state"`.
             * @details The system state will be saved to the file name `state_filename + "_%d.json"` at the specified
             * `save_frequency` as well as the first and last time increment.
             */
            std::string state_filename;

            std::string nodal_displacements_filename;
            std::string nodal_velocities_filename;
			std::string nodal_forces_filename;
        };

        ExplicitSystemManager(const std::string& config_filename);

		ExplicitSystem const& getExplicitSystem() const;

        /**
         *
         */
		rapidjson::Document const& getConfigDoc() const;

        /**
         * Returns time step used when integrating the equations of motion.
         */
		double getTimeStep() const;

        /**
         * @brief Returns the current number of times the system has been integrated.
         */
		unsigned int getIterationNumber() const;

        /**
         * @brief Runs the simulation
         * @details The equations of motion of the constructed system are integrated until the system time exceeds the
         * user-specifed `end_time`.
         */
        void run();

    private:
        void init_time_from_config();/**<Loads the start and end time from the configuration document.*/
        void construct_system();/**<Constructs an `ExplicitSystem` from the configuration document and estimates a stable time step.*/
        void dump_system();/**<Saves the current state of the system to a set of CSV files and JSON document.*/

        /**
         * @brief Saves a `ColumnVector` to the indicated file
         * @param vec Column vector to save to disk.
         * @param filename Name of file to save data to.
         */
        void save_column_vector(const ColumnVector& vec, std::string filename);

        Options options;/**<Set of user-specified options to customize the simulation.*/
        rapidjson::Document config_doc;/**<Parsed JSON document specified by the user. Parameters of the simulation.*/
        std::unique_ptr<ExplicitSystem> explicit_system;/**<System to be integrated. It is responsible for explicitly integration the equations of motion.*/
        double start_time;/**<The initial time of the system, i.e. the time that matches the time the user-specified initial conditions are valid.*/
        double end_time;/**<Simulation will stop when the current system time exceeds this value.*/
        double dt;/**<Time step. Amount of time that the system is advanced each iteration.*/
        unsigned int iteration_number;/**<Current iteration. The number of time steps the system has been integrated.*/
    };

} // namespace explicit_fea

#endif //EXPLICIT_FEA_EXPLICIT_SYSTEM_MANAGER_H
