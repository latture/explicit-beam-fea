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

#ifndef EXPLICIT_FEA_MESH_H
#define EXPLICIT_FEA_MESH_H

#include "beam_element.h"
#include <memory>
#include <set>
#include "bc.h"
#include "containers.h"

namespace explicit_fea {
    /**
     * A mesh is defined with a set of nodes, elements, and boundary conditions.
     * These inputs are used to construct global stiffness and mass matrices as
     * well as store the list of boundary conditions for use by the `ExplicitSystem`.
     */
    class Mesh {
    public:
        /**
         * @brief Default constructor
         * @details Creates an empty mesh
         */
        Mesh() {};

        /**
         * @brief Constructor
         * @details Builds the global stiffness and mass matrices as well as acting
         * as a source to store the boundary conditions as a member variable.
         * 
         * @param[in] nodes Node list containing a vector of \f$ (x, y, z) \f$ coordinates.
         * @param[in] elem List of elements referencing the nodal indices of the node list.
         * @param[in] bcs List of boundary conditions to apply to nodal degrees of freedom.
         */
        Mesh(const std::vector<Node>& nodes,
             const std::vector<std::unique_ptr<BeamElement>>& elems,
             BCList bcs);

        /**
         * @brief Returns the inverse mass matrix
         * @details The inverse mass matrix is returned as a sparse square matrix 
         * where the number of rows and columns equals the total number of nodal
         * degrees of freedom for the system.
         * 
         * @return Global inverse mass matrix
         */
        SparseMatrix const& get_inv_mass_matrix() const;

        SparseMatrix const& get_mass_matrix() const;

        /**
         * @brief Returns the global stiffness matrix
         * @details The stiffness matrix is returned as a sparse square matrix 
         * where the number of rows and columns equals the total number of nodal
         * degrees of freedom for the system.
         * 
         * @return Global stiffness matrix
         */        
        SparseMatrix const& get_global_stiffness_matrix() const;

        /**
         * @brief Returns the boundary conditions applied to the nodal degrees of freedom.
         * @return Vector of boundary conditions
         * @sa `BC`
         */
        BCList const& get_bcs() const;

    private:
        /**
         * @brief Assembles the global mass and stiffness matrices.
         * 
         * @param[in] nodes Node list containing a vector of \f$ (x, y, z) \f$ coordinates.
         * @param[in] elems Vector of elements referencing the associated node list.
         */
        void assemble_matrices(const std::vector<Node>& nodes, const std::vector<std::unique_ptr<BeamElement>>& elems);

        /**
         * @brief Converts the local elemental matrix into a set of triplets to later be used to construct a sparse global matrix.
         * 
         * @param[in] elem Current element
         * @param[in] triplets Set of triplets to append the elemental matrix data to.
         * @param[in] mat Local matrix associated with the current element.
         */
        void append_triplets(const BeamElement& elem,
                             std::vector<Eigen::Triplet<double> > &triplets,
                             const LocalMatrix& mat);
        
        /**
         * @brief Applies the boundary conditions to the global mass matrix.
         * @details The row and column of the inverse mass matrix associated
         * with a specified degree of freedom are affected by the boundary condition.
         * This method updates the assembled mass matrix to reflect the prescribed
         * nodal degrees of freedom.
         */
        void apply_bcs();

        void calc_mass_matrix();

        BCList bcs;/**<List of boundary conditions applied to nodal degrees of freedom*/
        SparseMatrix global_stiffness_matrix;/**<Assembled global stiffness matrix*/
        SparseMatrix mass_matrix;/**<Global mass matrix*/
        SparseMatrix inverse_mass_matrix;/**<Inverse of the global mass matrix*/

        /**
         * @brief Object used to prune the inverse mass matrix based on the list of boundary conditions.
         * @details This object facilitates straightforward removal of non-zero entries in the row and 
         * column associated with a prescribed degree of freedom.
         */
        struct BcsPruneFunctor {
            /**
             * @brief Constructor
             * @details Upon construction the list of boundary conditions are parsed into a set 
             * allowing efficient lookup when traversing the global inverse mass matrix in order
             * to determine if the current non-zero entries is affected by a boundary condition.
             * If the current entry is affected it is removed from the global sparse matrix.
             * 
             * @param[in] bcs List of boundary conditions applied to nodal degrees of freedom
             */
            BcsPruneFunctor(const BCList& bcs);

            /**
             * @brief Determines whether the current non-zero entry is affected by a boundary condition
             * and, thus, should be removed from the global mass matrix. 
             * 
             * @param[in] row Row of current non-zero entry.
             * @param[in] col Column of the current non-zero entry.
             * @param[in] value Value at the position specified by `inv_mass_matrix(row, col)`
             * @return Returns true if the value should be kept (i.e. not involved with a boundary condition) 
             * and false otherwise.
             */
            bool operator() (int row, int col, double value) const;

            std::set<int> bcs_ind_set;/**<Set of global degrees of freedom set by boundary conditions.*/
        };
    };
}

#endif //EXPLICIT_FEA_MESH_H
