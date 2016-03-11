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
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: ryan.latture@gmail.com (Ryan Latture)

#ifndef EXPLICIT_FEA_CONTAINERS_H
#define EXPLICIT_FEA_CONTAINERS_H

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <unsupported/Eigen/AlignedVector3>
#include <vector>

namespace explicit_fea {
    /**
     * An elemental matrix in local coordinates. Will either be the elemental stiffness matrix or the global-to-local rotation matrix
     */
    typedef Eigen::Matrix<double, 12, 12, Eigen::RowMajor> LocalMatrix;

    /**
     * Sparse matrix that is used internally to hold sparse representations of the global and elemental stiffness matrices.
     */
    typedef Eigen::SparseMatrix<double> SparseMatrix;

    /**
     * Vector that stores the nodal forces, i.e. the variable \f$[F]\f$ in \f$[K][Q]=[F]\f$,
     * where \f$[K]\f$ is the stiffness matrix and \f$[Q]\f$ contains the nodal displacements,
     * or nodal displacements  \f$[Q]\f$.
     */
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> ColumnVector;

    /**
     * @brief A node that describes a mesh. Uses Eigen's predefined Vector class for added functionality.
     * @details See the Eigen documentation on the `Vector3d` class for more options of what can be done with `Node`s. \n
     * Examples of constucting a `Node` at \f$(x, y, z)=(0,1,2)\f$:
     * @code
     * // specify values on constuction
     * explicit_fea::Node n1(1.0, 2.0, 3.0);
     *
     * // construct a Node then insert values
     * explicit_fea::Node n2;
     * n2 << 0.0, 1.0, 2.0;
     * @endcode
     */
    typedef Eigen::Vector3d Node;

    /**
     * @brief The set of properties associated with an element.
     * @details These properties are used with the elemental shape functions to fully define the 
     * elemental stiffness and mass matrices in an element formulation. These values are parsed,
     * for example, within the `EulerBernoulliBeamElement` to construct the extensional stiffness, \f$EA\f$, 
     * bending stiffness parallel to the local z-axis \f$EI_{z}\f$, bending stiffness parallel to the 
     * local y-axis\f$EI_{y}\f$, the torsional stiffness, \f$GJ\f$, among others.
     *
     * @code
     * double youngs_modulus = 200.0e9;
     * double shear_modulus = 80.0e9;
     * double area = 0.0314;
     * double Iz = 0.0002;
     * double Iy = 0.0002;
     * double J = 0.0004;
     * double density = 7800.0;
     * Eigen::Vector3d normal_vec;
     * normal_vec << 0.0, 0.0, 1.0;
     * fea::Props props(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, normal_vec);
     * @endcode
     */
    struct Props {
        const double youngs_modulus;/**<Young's (elastic) modulus.*/
        const double shear_modulus;/**<Shear modulus.*/
        const double area;/**<Cross-sectional area.*/
        const double Iz;/**<Second moment of area parallel to local z-axis.*/
        const double Iy;/**<Second moment of area parallel to local y-axis.*/
        const double J;/**<Torsional constant.*/
        const double density;/**<Density of the parent material.*/
        const Eigen::Vector3d normal_vec;/**<Vector normal to element (`size==3`). Direction should be parallel to the beam element's local y-axis.*/

        /**
         * @brief Constructor
         * @details Allows properties to be set upon initialization.
         *
         * @param[in] youngs_modulus Young's (elastic) modulus of parent material.
         * @param[in] shear_modulus Shear modulus of parent material.
         * @param[in] Iz Second moment of area parallel to local z-axis
         * @param[in] Iy Second moment of area parallel to local y-axis.
         * @param[in] J Torsional constant.
         * @param[in] density Density of the parent material.
         * @param[in] normal_vec Vector normal to element (`normal_vec.size()==3`). Direction should be parallel to the beam element's local y-axis.
         */
        Props(double _youngs_modulus, double _shear_modulus, double _area,
              double _Iz, double _Iy, double _J, double _density, const Eigen::Vector3d &_normal_vec)
                : youngs_modulus(_youngs_modulus), shear_modulus(_shear_modulus), area(_area),
                  Iz(_Iz), Iy(_Iy), J(_J), density(_density), normal_vec(_normal_vec) {
        };
    };

    /**
     * @brief Convenience enumerator for specifying the active degree of freedom in a constraint.
     */
    enum DOF {
        /**
         * Displacement along the global x-axis.
         */
                DISPLACEMENT_X,

        /**
         * Displacement along the global y-axis.
         */
                DISPLACEMENT_Y,

        /**
         * Displacement along the global z-axis.
         */
                DISPLACEMENT_Z,

        /**
         * Rotation about the global x-axis.
         */
                ROTATION_X,

        /**
         * Rotation about the global y-axis.
         */
                ROTATION_Y,

        /**
         * Rotation about the global z-axis.
         */
                ROTATION_Z,
        /**
         * Number of degrees of freedom per node.
         */
                NUM_DOFS
    };

} // namespace fea

#endif // EXPLICIT_FEA_CONTAINERS_H
