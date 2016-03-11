// Copyright 2016. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, ORjavascript:void(0)
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: ryan.latture@gmail.com (Ryan Latture)

#ifndef EXPLICIT_FEA_SYSTEM_H
#define EXPLICIT_FEA_SYSTEM_H

#include "mesh.h"
#include <Eigen/SparseCholesky>
#include <memory>
#include <rapidjson/document.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/stringbuffer.h>
#include <string>
#include "containers.h"
#include "fe_utils.h"
#include "force.h"

namespace explicit_fea {
    /**
     * @brief System responsible managing explicit time integration.
     * @details Each system holds a mesh and set of external forces
     * and advancing the solution in time by explicitly integrating
     * the equations of motion using the Newmark-&beta; method.
     */
    class ExplicitSystem {
    public:

        /**
         * @brief Options used to customize the behavior of the system.
         */
        struct Options {
            Options() :
                    beta(0.25),
                    gamma(0.5),
                    damping_beta(0.01),
                    damping_alpha(0.01)
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
             * &beta; constant referred to in the Newmark-&beta; method, Default = `0.25`.
             */
            double beta;

            /**
             * &gamma; constant referred to in the Newmark-&beta; method, Default = `0.5`.
             */
            double gamma;

            /**
             * Damping constant &alpha; used to construct the proportional damping matrix, Default = `0.01`.
             * The damping matrix is formulated from the equation:
             * \f$ \left[C \right] = \alpha \left[ M \right] + \beta \left[ K \right] \f$ \n
             * where \f$ \left[ M \right] \f$ is the global mass matrix,
             * \f$ \left[ K \right] \f$ is the global stiffness matrix,
             * and \f$ \left[ C \right] \f$ is the proportional damping matrix.
             */
            double damping_alpha;

            /**
             * Damping constant &beta; used to construct the proportional damping matrix, Default = `0.01`.
             * The damping matrix is formulated from the equation:
             * \f$ \left[C \right] = \alpha \left[ M \right] + \beta \left[ K \right] \f$ \n
             * where \f$ \left[ M \right] \f$ is the global mass matrix,
             * \f$ \left[ K \right] \f$ is the global stiffness matrix,
             * and \f$ \left[ C \right] \f$ is the proportional damping matrix.
             */
            double damping_beta;
        };

        /**
         * @brief Default constructor
         * @details Initializes an empty system.
         */
        ExplicitSystem() {};

        /**
         * @brief Constructor
         * @details Constructs a system containing the specified mesh and applied 
         * external forces. To avoid copy construction of what, in general, are 
         * expensive objects both the mesh and set of external forces must be 
         * cast as r-values to be "moved" into the system. In effect the explicit
         * system acts as a sink for both input parameters.
         *
         * @param[in] mesh. Mesh constructed from a set of node and element lists
         *              in coordination with a set of boundary conditions. Responsible for
         *              providing the global stiffness and mass matrices.
         * @param[in] external_forces. Vector of nodal forces applied to the model.
         * @param[in] initial_displacements Initial nodal displacements.
         * @param[in] initial_velocities Initial nodal velocities.
         * @param[in] t0 Initial time of the system.
         * @param[in] options Options used to customize behavior of the system.
         */
        ExplicitSystem(Mesh mesh,
                       ForceList external_forces,
                       const ColumnVector& initial_displacements,
                       const ColumnVector& initial_velocities,
                       double t_0,
                       const Options& options);

        /**
         * @brief Updates the system to the specified time.
         * @details Explicit time integration is used to advance the system from the
         * last previous time to that specified by the input parameter.
         * 
         * @param[in] dt Time increment. Amount of time to advance the system by.
         */
        void update(double dt);

        /**
         * @brief Returns the current nodal displacements.
         * @details Nodal displacements are held in a 1D vector of the form \n
         * \f$ \left[ 
         *     \delta u_{1}^{x}, \delta u_{1}^{y}, \delta u_{1}^{z},
         *     \delta \theta_{1}^{x}, \delta \theta_{1}^{y}, \delta \theta_{1}^{z},
         *     ...,
         *     \delta u_{N}^{x}, \delta u_{N}^{y}, \delta u_{N}^{z},
         *     \delta \theta_{N}^{x}, \delta \theta_{N}^{y}, \delta \theta_{N}^{z},
         *   \right] \f$ \n
         * where \f$ \delta u \f$ specifies displacement along translational degrees
         * of freedom for the given node number \f$ \left( 1, 2, 3, ...N\right) \f$ and 
         * direction \f$ \left( x, y, z) \f$ and \f$ \delta \theta \f$ specifies the
         * change in rotation (in radians) of the noted nodal degree of freedom.
         *   
         * @return Nodal displacements at the current time.
         */
        ColumnVector const& getDisplacements() const;

        /**
         * @brief Returns the current nodal forces.
         * @details Nodal forces are held in a 1D vector of the form \n
         * \f$ \left[ 
         *     P_{1}^{x}, P_{1}^{y}, P_{1}^{z},
         *     M_{1}^{x}, M_{1}^{y}, M_{1}^{z},
         *     ...,
         *     P_{N}^{x}, P_{N}^{y}, P_{N}^{z},
         *     M_{N}^{x}, M_{N}^{y}, M_{N}^{z},
         *   \right] \f$ \n
         * where \f$ P \f$ specifies nodal force along translational degrees
         * of freedom for the given node number \f$ \left( 1, 2, 3, ...N\right) \f$ and 
         * direction \f$ \left( x, y, z) \f$ and \f$ \delta \theta \f$ specifies the
         * moment of the noted nodal degree of freedom.
         *   
         * @return Nodal forces at the current time.
         */
        ColumnVector getForces() const;

        /**
         * @brief Returns the current nodal velocities.
         * @details Nodal velocities are held in a 1D vector of the form \n
         * \f$ \left[ 
         *     v_{1}^{x}, v_{1}^{y}, v_{1}^{z},
         *     \omega_{1}^{x}, \omega_{1}^{y}, \omega_{1}^{z},
         *     ...,
         *     v_{N}^{x}, v_{N}^{y}, v_{N}^{z},
         *     \omega_{N}^{x}, \omega_{N}^{y}, \omega_{N}^{z},
         *   \right] \f$ \n
         * where \f$ v \f$ specifies nodal velocity along translational degrees
         * of freedom for the given node number \f$ \left( 1, 2, 3, ...N\right) \f$ and 
         * direction \f$ \left( x, y, z) \f$ and \f$ \delta \theta \f$ specifies the
         * angular velocity of the noted nodal degree of freedom.
         *   
         * @return Nodal velocities at the current time.
         */
        ColumnVector const& getVelocities() const;

        /**
         * @brief Returns the current time of the system.
         */
        double getTime() const;

        /**
         * Returns the `Mesh` contained within the system.
         */
		Mesh const& getMesh() const;

    private:

        /**
         * @brief Creates the proportional damping matrix.
         * @details The proportional damping matrix is given by:
         * \f$ \alpha \left[ M \right] + \beta \left[ K \right] \f$, where \f$ \alpha \f$
         * and \f$ \beta \f$ are constants specified by in the
         * `ExplicitSystem::Options` struct, \f$ \left[ M \right]  \f$ is the global mass matrix,
         * and \f$ \left[ K \right]  \f$ is the global stiffness matrix.
         */
        void assemble_damping_matrix();

        /**
         * @brief Assembles the left hand side of the equation of motion
         * @details The left-hand side of the equation of motion is given by: \n
         * \f$ \left[ M \right] + \gamma \Delta t \left[ C \right] + \beta \Delta t^{2} \left[ K \right] \f$ \n
         * and needs to be reassembled every time the current time step does not match the previous one.
         */
        void assemble_LHS(double dt);

        /**
         * Applies the boundary conditions of the mesh to the nodal
         * displacements, velocities, and accelerations.
         * @param[in] time Time to advance the system to.
         */
        void apply_bcs(double time);

        /**
         * Applies the external forces applied to nodal degrees of freedom.
         * @param[in] time Time to advance the system to.
         */
        void apply_external_forces(double time);

        Mesh mesh;/**<Mesh constructed from a set of node and element lists
                    * in coordination with a set of boundary conditions. Responsible for
                    * providing the global stiffness and mass matrices.*/
        ForceList external_forces;/**<Time-variant external nodal forces applied to the mesh.*/
        ColumnVector displacements_0;/**<Current nodal displacements.*/
        ColumnVector velocities_0;/**<Current nodal velocities.*/
        ColumnVector velocities_1;/**<Nodal velocities for the next iteration's (\f$ i+1 \f$) system time.*/
        ColumnVector accelerations_0;/**<Current nodal accelerations.*/
        ColumnVector accelerations_1;/**<Nodal accelerations for the next iteration's (\f$ i+1 \f$) system time.*/
        ColumnVector forces;/**<Current nodal forces.*/
        ColumnVector RHS;/**<Effective nodal forces.*/
        double t_0;/**<Current time of the system.*/
        double dt_0;/**Time step used for last integration update.*/
        const Options options;/**<User specified options to modify the behavior of the system.*/
        SparseMatrix damping_matrix;/**<Proportional damping matrix.*/
        /**
         * Left hand side of the equation of motion given by: \n
         * \f$ \left[ M \right] + \gamma \Delta t \left[ C \right] + \beta \Delta t^{2} \left[ K \right] \f$. \n
         * This is assembled every time the supplied time step does not match the previous \f$ \Delta t \f$.
         */
        SparseMatrix LHS;
        Eigen::SimplicialLDLT<SparseMatrix> solver;/**<Solver used to calculate the accelerations at the next time based off */
        ValueCompare<double> compare;/**<Object used to compare the time steps in order to reassemble the left hand side of the equations of motion.*/
    };

} // namespace explicit_fea

#endif // EXPLICIT_FEA_SYSTEM_H