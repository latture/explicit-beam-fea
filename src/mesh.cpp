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

#include "mesh.h"
#include <Eigen/Dense>
#include <iostream>

namespace explicit_fea {
    Mesh::Mesh(const std::vector<Node> &nodes,
               const std::vector<std::unique_ptr<BeamElement>>& elems,
               BCList bcs) :
        bcs(std::move(bcs))
    {
        assemble_matrices(nodes, elems);
        apply_bcs();
    }
    
    SparseMatrix const& Mesh::get_global_stiffness_matrix() const {
        return global_stiffness_matrix;
    }
    
    SparseMatrix const& Mesh::get_inv_mass_matrix() const {
        return inverse_mass_matrix;
    }

    SparseMatrix const& Mesh::get_mass_matrix() const {
        return mass_matrix;
    }

    BCList const& Mesh::get_bcs() const {
        return bcs;
    }

    void Mesh::assemble_matrices(const std::vector<Node> &nodes, const std::vector<std::unique_ptr<BeamElement>>& elems) {
        int n = DOF::NUM_DOFS * nodes.size();
        global_stiffness_matrix.resize(n, n);
        mass_matrix.resize(n, n);
        inverse_mass_matrix.resize(n, n);

        std::vector<Eigen::Triplet<double> > stiffness_triplets, mass_triplets, inv_mass_triplets;
        stiffness_triplets.reserve(40 * elems.size());
        mass_triplets.reserve(40 * elems.size());
        inv_mass_triplets.reserve(40 * elems.size());
        LocalMatrix inv_mass, mass;

        for (size_t i = 0; i < elems.size(); ++i) {
            inv_mass = elems[i]->calculate_inv_mass_matrix(nodes);
            mass = inv_mass.inverse();
            append_triplets(*elems[i], stiffness_triplets, elems[i]->calculate_stiffness_matrix(nodes));
            append_triplets(*elems[i], mass_triplets, mass);
            append_triplets(*elems[i], inv_mass_triplets, inv_mass);
        }
        global_stiffness_matrix.setFromTriplets(stiffness_triplets.begin(), stiffness_triplets.end());
        inverse_mass_matrix.setFromTriplets(inv_mass_triplets.begin(), inv_mass_triplets.end());
        mass_matrix.setFromTriplets(mass_triplets.begin(), mass_triplets.end());

        global_stiffness_matrix.prune(1.e-14);
        mass_matrix.prune(1.e-14);
        inverse_mass_matrix.prune(1.e-14);

        global_stiffness_matrix.makeCompressed();
        mass_matrix.makeCompressed();
        inverse_mass_matrix.makeCompressed();
    }

    void Mesh::append_triplets(const BeamElement& elem, std::vector<Eigen::Triplet<double> > &triplets, const LocalMatrix& mat) {
        const int nn1 = elem.get_node_numbers()[0];
        const int nn2 = elem.get_node_numbers()[1];
        size_t row, col;
        const unsigned int dofs_per_elem = DOF::NUM_DOFS;

        SparseMatrix sparse_mat = mat.sparseView();

        for (size_t j = 0; j < sparse_mat.outerSize(); ++j) {
            for (SparseMatrix::InnerIterator it(sparse_mat, j); it; ++it) {
                row = it.row();
                col = it.col();

                // check position in local matrix and update corresponding global position
                if (row < dofs_per_elem) {
                    // top left
                    if (col < dofs_per_elem) {
                        triplets.push_back(Eigen::Triplet<double>(dofs_per_elem * nn1 + row,
                                                                  dofs_per_elem * nn1 + col,
                                                                  it.value()));
                    }
                        // top right
                    else {
                        triplets.push_back(Eigen::Triplet<double>(dofs_per_elem * nn1 + row,
                                                                  dofs_per_elem * (nn2 - 1) + col,
                                                                  it.value()));
                    }
                }
                else {
                    // bottom left
                    if (col < dofs_per_elem) {
                        triplets.push_back(Eigen::Triplet<double>(dofs_per_elem * (nn2 - 1) + row,
                                                                  dofs_per_elem * nn1 + col,
                                                                  it.value()));
                    }
                        // bottom right
                    else {
                        triplets.push_back(Eigen::Triplet<double>(dofs_per_elem * (nn2 - 1) + row,
                                                                  dofs_per_elem * (nn2 - 1) + col,
                                                                  it.value()));
                    }
                }
            }
        }
    }

    void Mesh::apply_bcs() {
        BcsPruneFunctor prune_functor(bcs);
        inverse_mass_matrix.prune(prune_functor);
        for (std::set<int>::const_iterator it = prune_functor.bcs_ind_set.begin();
             it != prune_functor.bcs_ind_set.end(); ++it) {
            inverse_mass_matrix.insert(*it, *it) = 1.0;
        }
    }

    Mesh::BcsPruneFunctor::BcsPruneFunctor(const BCList &bcs){
        for (BCList::const_iterator it = bcs.begin(); it != bcs.end(); ++it) {
            bcs_ind_set.insert((*it)->global_index);
        }
    }

    bool Mesh::BcsPruneFunctor::operator()(int row, int col, double value) const {
        bool row_exists = bcs_ind_set.find(row) != bcs_ind_set.end();
        if (row_exists)
            return false;

        bool col_exists = bcs_ind_set.find(col) != bcs_ind_set.end();
        return !col_exists;
    }

} // namespace explicit_fea