# 3D linear beam element code  #

## Getting started ##
This project requires CMake to compile the code. If not installed, please install CMake before continuing.

#### To compile the code: ####
  1. Open the `explicit-beam-fea` directory in the terminal
  2. Execute `git submodule update --init --recursive` to fetch all the dependencies.
  2. Create a folder named `build` within the `explicit-beam-fea` parent directory.
  3. Open a terminal and navigate to the newly formed `build` directory
  4. Execute `cmake ..`
    * Use `-DCMAKE_BUILD_TYPE=debug` if you would like to build the code for debugging purposes. By default the make files will be configured for the release build.
  5. On Linux run `make` in the terminal from the build directory to build all the targets. On Windows the solution file will be located in the build directory. Open the solution file in Visual Studio and compile.

## Introduction ##
This contains a dynamic finite element code intended for either 3D Euler-Bernoulli or Timoshenko beam elements.
Advancing the system is accomplished via explicit time integration of the governing equations.
Formulating an analysis can be done in C++, or if time-independent boundary conditions are sufficient the command 
line interface can be used.

### Method 1: Using C++ ###

An analysis consists a node list (`explicit_fea::Node`), element list (`explicit_fea::BeamElement`), boundary conditions (`explicit_fea::BC`), 
prescribed nodal forces (`explicit_fea::Force`), and initial conditions (specifically nodal velocities and displacements).
The nodes, elements, and boundary conditions are parsed into a mesh (`explicit_fea::Mesh`), which is in turn responsible for assembling the global stiffness and mass matrices.
The mesh is added to an explicit system which is responsible for advancing the solution via explicit time incrementation.

#### Forming the mesh ####
The mesh defines the nodal coordinates in `(x, y, z)` space, nodes that are connected to form beam elements, elemental properties, 
and boundary conditions applied to the mesh (as they will affect the formulation of the mass matrix).
The nodal coordinates are formed as a vector of `explicit_fea::Node`'s where each node simply contains the `(x, y, z)` coordinates of the point.
An element contains the 2 nodal indices that are connected to form the element as well as the associated properties of the element.
The properties must define the Young's modulus of the parent material (`youngs_modulus`), shear modulus of the parent material (`shear_modulus`), 
cross-sectional area (`area`), second moment of area parallel to the local z-axis `EIz`,
second moment of area parallel to the local y-axis (`EIy`), torsion constant (`J`), and a vector pointing along the beam elements local y-axis.
Boundary conditions specify prescribed nodal velocities or displacements. Each must specify the node number, degree of freedom, type, and 
implement a function that returns the prescribed value at a given time. The latter requires that boundary conditions inherit from the `explicit_fea::BC`
class and implement the `double get_value(double time)` function. An example of a boundary condition that is constant with respect to time has been
implemented as the `explicit_fea::ConstantBC` class. The index of the node is simply the index the node occurs in the node list. 
The degree of freedom can be defined using the `explicit_fea::DOF` enum or by specifying the integer associated with the degree of freedom explicitly. 
There are 6 degrees of freedom per node meaning valid integers associated with degrees of freedom are between 0 and 5. 
The associations for degrees of freedom are defined as

  * 0 = DISPLACEMENT_X, displacement along the global x-axis.
  * 1 = DISPLACEMENT_Y, displacement along the global y-axis.
  * 2 = DISPLACEMENT_Z, displacement along the global z-axis.
  * 3 = ROTATION_X, rotation about the global x-axis (in radians).
  * 4 = ROTATION_Y, rotation about the global y-axis (in radians).
  * 5 = ROTATION_Z, rotation about the global z-axis (in radians).

Valid boundary condition types are specified in `explicit_fea::BC::Type` enum. There are currently two options

  * 0 = DISPLACEMENT
  * 1 = VELOCITY

An example forming a simple mesh with a single element is shown below.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
// [form the node list
explicit_fea::Node node1, node2;

// place the first node at (0, 0, 0)
node1 << 0, 0, 0;

// place the second node at (1, 0, 0)
node2 << 1, 0, 0;

std::vector<explicit_fea::Node> node_list = {node1, node2};
// ]

// [ form the element list
// define the indices of the node list that form the element
unsigned int nn1 = 0;
unsigned int nn2 = 1;

// define the properties of the element
double youngs_modulus = 200.0e9
double shear_modulus = 80.0e9
double area = 0.112;
double Iz = 0.001;
double Iy = 0.001;
double J = 0.002;
Eigen::Vector3d normal_vec;
normal_vec << 0.0, 0.0, 1.0;
fea::Props props(youngs_modulus, shear_modulus, area
                 Iz, Iy, J, normal_vec);

// The element list must be a vector of pointers to a
// beam element implementation. The currently supported
// beam elements are TimoshenkoBeamElement and 
// EulerBernoulliBeamElement
std::vector<std::unique_ptr<explicit_fea::Elem>> elem_list;
elem_list.push_back(std::unique_ptr<TimoshenkoBeamElement>(new TimoshenkoBeamElement(nn1, nn2, props)));
// ]

// [ define boundary conditions

// a BCList is simply a typedef for std::vector<std::unique_ptr<BC>>.
BCList bcs;

// fix all of node1's displacement degrees of freedom
for (unsigned int i = 0; i < explicit_fea::DOF::NUM_DOFS, ++i)
    bcs.push_back(std::unique_ptr<ConstantBC>(new ConstantBC(0, i, 0.0, explicit_fea::BC::Type::DISPLACEMENT)))

// Add a positive velocty on the x degree of freedom of node2
bcs.push_back(std::unique_ptr<ConstantBC>(new ConstantBC(1, 0, 0.01, explicit_fea::BC::Type::VELOCITY)))
// ]

// [ construct the mesh
// nodes and elems are passed by const reference
// boundary conditions must be moved into the mesh
// because they are held as a member variable and 
// are needed by the ExplicitSystem later.
explicit_fea::Mesh mesh(nodes, elems, std::move(bcs))
// ]
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Nodal forces ####
Nodal forces are external forces applied to nodal degrees of freedom.
They are assigned in a similar manner to boundary conditions but there is no type specification, 
i.e. only node number, degree of freedom, and value. Prescribed nodal forces must inherit from
`explicit_fea::Force` and implement the `double get_value(double time)` method, which returns
the value of the degree of freedom at a given time. The `explicit_fea::ConstantForce` class
has been defined for external forces that do not vary with time.
We can load our previously defined mesh in the global y-direction using:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
// define the node index and value
unsigned int nn2 = 1;
double value = 100.0;

// A ForceList is simply a typedef for std::vector<std::unique_ptr<Force>>
ForceList force_list;

// create the force and add to the vector
force_list.push_back(std::unique_ptr<ConstantForce>(new ConstantForce(nn2, explicit_fea::DOF::DISPLACEMENT_Y, value)))
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Explicit system ####
The explicit system holds the mesh and external forces and is responsible for integrating the equations of motion.
It accomplishes this using the Newmark-beta method. In addition to the mesh and forces, the system requires the
initial conditions to be specified. These are the initial nodal displacements, nodal velocities, and system time.
An example of constructing an explicit system is shown below

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
// see the previous section on how to construct the mesh and forces.
// let's assume that we have previously constructed these values
// for our single element system.
explicit_fea::Mesh mesh;
ForceList forces;

// we can customize the behavior of the system by using the options struct
explicit_fea::ExplicitSystem::Options options; // creates an object with the default behavior
// we can change the beta and gamma coefficents to use the linear acceleration method:
options.beta = ;
options.gamma = ;
// let's also change the parameters for the proportional damping matrix
options.damping_alpha = 0.001;
options.damping_beta = 0.001;

// [ define the initial conditions
// each column vector should have the same number of rows as the global matrices
auto num_cols = mesh.get_stiffness_matrix().cols();

// initialize displacements and velocities with values of 0 for all entries
explicit_fea::ColumnVector initial_displacements(explicit_fea::ColumnVector::Zero(num_cols));
explicit_fea::ColumnVector initial_velocities(explicit_fea::ColumnVector::Zero(num_cols));

// time corresponding the initial displacements and velocities
double initial_time = 0.0;
// ]

// now we can construct the system:
// we move the mesh and forces into the system as 
// they are stored as member variables and contain
// unique_ptr's and cannot be copy constructed
explicit_fea::ExplicitSystem es(std::move(mesh), std::move(forces), 
                                initial_displacements, initial_velocities,
                                initial_time, options)

// advancing the system state can be done using:
double dt = 1.e-6;
es.update(dt)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#### Managing a simulation ####

The `ExplicitSystemManager` was created to help instantiate and run a simulation.
It parses a configuration document, constructs a system, estimates a stable time
step, then integrates the system until the desired end time is exceeded. Using 
this class is appropriate when external forces and boundary conditions are constant
in time. Constructing and running a simulation resembles

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
// define the filename containing the configuration file
std::string config_file = "/path/to/config.json"
explicit_fea::ExplicitSystemManager manager(config_filename)
manager.run()
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The configuration file is a json document that has key value pair associated with 
the information needed to construct the system and itialize state. An example is 
listed below

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
{
    "nodes"   : "path/to/nodes.csv",
    "elems"   : "path/to/elems.csv",
    "props"   : "path/to/props.csv",
    "bcs"     : "path/to/bcs.csv",
    "forces"  : "path/to/forces.csv",
    "nodal_displacements"  : "path/to/nodal_displacements.csv",
    "nodal_velocities"     : "path/to/nodal_velocities.csv",
    "start_time" : 0.0,
    "end_time"   : 1.0,
    "options" : {
                    "verbose" : true,
                    "save_frequency" : 10,
                    "state_filename" : "state",
                    "nodal_displacements_filename" : "displacements",
                    "nodal_velocities_filename" : "velocities",
                    "nodal_forces_filename" : "forces",
                    "beta" : 0.25,
                    "gamma" : 0.5,
                    "damping_alpha" : 0.01,
                    "damping_beta" : 0.01
                }
}
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Formatting the CSV files is described in the following section.
Each key must have a specific format in order to be parsed.
For more information on the optional parameters see the 
`explicit_fea::ExplicitSystemManager::Options` and `explicit_fea::ExplicitSystem::Options`
classes. Options for both are given in the same key in the configuration file and parsed
as-needed when constructing the system and system manager.

The use of a JSON document avoids the need to set each of these options using command line options, which can become tedious when running multiple jobs.
The "nodes", "elems", "props", "nodal_displacements", "nodal_velocites", "start_time", and "end_time" keys are required. 
Keys "bcs" and "forces" are optional--if not provided the analysis will assume none were prescribed.
If the "options" key is not provided the analysis will run with the default options.
Any of all of the "options" keys presented above can be used to customize the analysis.
If a key is not provided the default value is used in its place.
See the Formatting CSV Files section below for how the CSV files should be created.

#### Formatting CSV files ####
All CSV file must be comma delimited with no spaces between values, i.e. one row of the nodal coordinates file might resemble `1.0,2.0,3.0`.
The file indicated by the value of "nodes" should be in the format:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
x1,y1,z1
x2,y2,z2
...
...
...
xN,yN,zN
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

where each entry is a double and every line must have 3 entries for the `x,y,z` position.
The "elems" file contains (only) the node indices:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
el1_node_num_1,el1_node_num_2
el2_node_num_1,el2_node_num_2
...
...
...
elN_node_num_1,elN_node_num_2
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

where each entry is an integer and must have 2 nodal indices defining the connectivity of the element.
Elemental properties are defined in the "props" file as:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
el1_youngs_modulus,el1_shear_modulus,el1_area,el1_Iz,el1_Iy,el1_J,el1_nvec_x_comp,el1_nvec_y_comp,el1_nvec_z_comp
el2_youngs_modulus,el2_shear_modulus,el2_area,el2_Iz,el2_Iy,el2_J,el2_nvec_x_comp,el2_nvec_y_comp,el2_nvec_z_comp
...
...
...
elN_youngs_modulus,elN_shear_modulus,elN_area,elN_Iz,elN_Iy,elN_J,elN_nvec_x_comp,elN_nvec_y_comp,elN_nvec_z_comp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

where each entry is a double and each line has 10 entries. The number of rows in the properties file must equal
the total number of elements in the mesh.

Each line in the "bcs" file specifies the node number, degree of freedom, value, and type:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
bc1_node_num,bc1_dof,bc1_value,bc1_type
bc2_node_num,bc2_dof,bc2_value,bc2_type
...
...
...
bcN_node_num,bcN_dof,bcN_value,bcN_type
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

where the node number is the index of the node in the node list (integer),
the DOF is the degree of freedom constrained (integer between 0 and 5),
value is the value to hold the degree of freedom at relative to the starting position (double),
and type is either 0 or 1. Type 0 indicates that the boundary condition applies to the displacement
degree of freedom and 1 indicates that is applies to velocity.

The "forces" CSV file is specified using the format:

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.txt}
force1_node_num,force1_dof,force1_value
force2_node_num,force2_dof,force2_value
...
...
...
forceN_node_num,forceN_dof,forceN_value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### Method 2: Using the command line interface ###
After using CMake to build the targets, an executable will be created that provide a command line interface (CLI) to the beam element code.
Once in the build directory, navigate to the `bin` folder containing fea_cmd. running `./explicit_fea_cmd -h` from the terminal will show the help
documentation for the program. The CLI expects the `-c` flag to be set and point to the config file for the current analysis.
A config file is a JSON document that contains key, value pairs pointing to the nodes, elements, properties, and other analysis options.
Internally this simply constructs an `ExplicitSystemManager` from the file and runs the analysis.

## Contact ##

* Ryan Latture (ryan.latture@gmail.com)
