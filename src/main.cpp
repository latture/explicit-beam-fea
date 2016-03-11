
#include <tclap/CmdLine.h>
#include "explicit_system_manager.h"

int main(int argc, char *argv[]) {
    try {
        TCLAP::CmdLine cmd("3D explicit beam element FEA. "
                                   "Use the -c [--config] flag to point to the configuration file for the current analysis.",
                           ' ', "1.0");
        TCLAP::ValueArg<std::string> configArg("c",
                                               "config",
                                               "Finite element configuration file (json format). "
                                                       "Must have \"nodes\", \"elems\", \"props\", \"nodal_displacements\", and "
                                                       "\"nodal_velocities\", members pointing the associated files. Keys "
                                                       "\"start_time\" and \"end_time\" must specify the time period of the analysis." 
                                                       "Optionally, any boundary conditions should be in the file pointed to by "
                                                       "the \"bcs\" member variable of the config file. Likewise, prescribed "
                                                       "forces are set using the \"forces\" variable. Please refer to the "
                                                       "documentation for the file format of each variable. Override the default "
                                                       "options using the \"options\" member variable the itself is a nested json object. "
                                                       "Refer to the documentation the possible configurations that can be set.",
                                               true,
                                               "config.json",
                                               "string");
        cmd.add(configArg);
        cmd.parse(argc, argv);
        std::string config_filename = configArg.getValue();
        explicit_fea::ExplicitSystemManager manager(config_filename);
        manager.run();
		std::cout << "Analysis completed." << std::endl;
    }
    catch (TCLAP::ArgException &e)  // catch any exceptions from parsing
    {
        std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
    }
    catch (std::exception &e) {
        std::cerr << "error: " << e.what() << std::endl;
    }
    return 0;
}