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

#include "fe_utils.h"
#include <Eigen/Dense>
#include <limits>
#include <memory>

namespace explicit_fea {

    void update_rotation(const std::vector<Node> &nodes, const BeamElement& elem, LocalMatrix &rotation, LocalMatrix &rotation_transposed) {
        rotation.setZero();
        rotation_transposed.setZero();
        const int nn1 = elem.get_node_numbers()[0];
        const int nn2 = elem.get_node_numbers()[1];

        // calculate unit normal vector along local x-direction
        Eigen::Vector3d nx = nodes[nn2] - nodes[nn1];
        nx.normalize();

        // calculate unit normal vector along y-direction
        const Eigen::Vector3d ny = elem.get_props().normal_vec.normalized();

        // calculate the unit normal vector in local z direction
        Eigen::Vector3d nz = nx.cross(ny);
        nz.normalize();

        // update rotation matrix
        rotation(0, 0) = nx(0);
        rotation(0, 1) = nx(1);
        rotation(0, 2) = nx(2);
        rotation(1, 0) = ny(0);
        rotation(1, 1) = ny(1);
        rotation(1, 2) = ny(2);
        rotation(2, 0) = nz(0);
        rotation(2, 1) = nz(1);
        rotation(2, 2) = nz(2);

        rotation(3, 3) = nx(0);
        rotation(3, 4) = nx(1);
        rotation(3, 5) = nx(2);
        rotation(4, 3) = ny(0);
        rotation(4, 4) = ny(1);
        rotation(4, 5) = ny(2);
        rotation(5, 3) = nz(0);
        rotation(5, 4) = nz(1);
        rotation(5, 5) = nz(2);

        rotation(6, 6) = nx(0);
        rotation(6, 7) = nx(1);
        rotation(6, 8) = nx(2);
        rotation(7, 6) = ny(0);
        rotation(7, 7) = ny(1);
        rotation(7, 8) = ny(2);
        rotation(8, 6) = nz(0);
        rotation(8, 7) = nz(1);
        rotation(8, 8) = nz(2);

        rotation(9, 9) = nx(0);
        rotation(9, 10) = nx(1);
        rotation(9, 11) = nx(2);
        rotation(10, 9) = ny(0);
        rotation(10, 10) = ny(1);
        rotation(10, 11) = ny(2);
        rotation(11, 9) = nz(0);
        rotation(11, 10) = nz(1);
        rotation(11, 11) = nz(2);

        // update transposed rotation matrix
        rotation_transposed(0, 0) = nx(0);
        rotation_transposed(0, 1) = ny(0);
        rotation_transposed(0, 2) = nz(0);
        rotation_transposed(1, 0) = nx(1);
        rotation_transposed(1, 1) = ny(1);
        rotation_transposed(1, 2) = nz(1);
        rotation_transposed(2, 0) = nx(2);
        rotation_transposed(2, 1) = ny(2);
        rotation_transposed(2, 2) = nz(2);

        rotation_transposed(3, 3) = nx(0);
        rotation_transposed(3, 4) = ny(0);
        rotation_transposed(3, 5) = nz(0);
        rotation_transposed(4, 3) = nx(1);
        rotation_transposed(4, 4) = ny(1);
        rotation_transposed(4, 5) = nz(1);
        rotation_transposed(5, 3) = nx(2);
        rotation_transposed(5, 4) = ny(2);
        rotation_transposed(5, 5) = nz(2);

        rotation_transposed(6, 6) = nx(0);
        rotation_transposed(6, 7) = ny(0);
        rotation_transposed(6, 8) = nz(0);
        rotation_transposed(7, 6) = nx(1);
        rotation_transposed(7, 7) = ny(1);
        rotation_transposed(7, 8) = nz(1);
        rotation_transposed(8, 6) = nx(2);
        rotation_transposed(8, 7) = ny(2);
        rotation_transposed(8, 8) = nz(2);

        rotation_transposed(9, 9) = nx(0);
        rotation_transposed(9, 10) = ny(0);
        rotation_transposed(9, 11) = nz(0);
        rotation_transposed(10, 9) = nx(1);
        rotation_transposed(10, 10) = ny(1);
        rotation_transposed(10, 11) = nz(1);
        rotation_transposed(11, 9) = nx(2);
        rotation_transposed(11, 10) = ny(2);
        rotation_transposed(11, 11) = nz(2);
    }

    double estimate_stable_timestep(const std::vector<Node>& nodes, const std::vector<std::unique_ptr<BeamElement>>& elems) {
        double length, timestep, wavespeed;
        double min_timestep = std::numeric_limits<double>::max();
        Eigen::Vector3d dn;
        int nn1, nn2;

        for (size_t i = 0; i < elems.size(); ++i) {
            nn1 = elems[i]->get_node_numbers()[0];
            nn2 = elems[i]->get_node_numbers()[1];

            dn = nodes[nn2] - nodes[nn1];
            length = dn.norm();
            wavespeed = sqrt(elems[i]->get_props().youngs_modulus / elems[i]->get_props().density);
            timestep = length / wavespeed;
            if (timestep < min_timestep) min_timestep = timestep;
        }
        return min_timestep / 10.0;
    }

} // namespace explicit_fea