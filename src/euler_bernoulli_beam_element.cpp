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

#include "euler_bernoulli_beam_element.h"
#include "fe_utils.h"

namespace explicit_fea {

    EulerBernoulliBeamElement::EulerBernoulliBeamElement(int nn1, int nn2, const Props& props) : BeamElement(nn1, nn2, props) {
    }

    LocalMatrix EulerBernoulliBeamElement::calculate_stiffness_matrix(const std::vector <Node> &nodes) {
        // extract element properties
        const double EA = props.youngs_modulus * props.area;
        const double EIz = props.youngs_modulus * props.Iz;
        const double EIy = props.youngs_modulus * props.Iy;
        const double GJ = props.shear_modulus * props.J;

        // store node indices of current element
        const int nn1 = node_numbers[0];
        const int nn2 = node_numbers[1];

        // calculate the length of the element
        const Node dn = nodes[nn1] - nodes[nn2];
        const double length = dn.norm();

        // store the entries in the (local) elemental stiffness matrix as temporary values to avoid recalculation
        const double tmpEA = EA / length;
        const double tmpGJ = GJ / length;

        const double tmp12z = 12.0 * EIz / (length * length * length);
        const double tmp6z = 6.0 * EIz / (length * length);
        const double tmp1z = EIz / length;

        const double tmp12y = 12.0 * EIy / (length * length * length);
        const double tmp6y = 6.0 * EIy / (length * length);
        const double tmp1y = EIy / length;

        LocalMatrix local_elemental_stiffness, rotation, rotation_transposed;
        local_elemental_stiffness.setZero();

        // update local elemental stiffness matrix
        local_elemental_stiffness(0, 0) = tmpEA;
        local_elemental_stiffness(0, 6) = -tmpEA;

        local_elemental_stiffness(1, 1) = tmp12z;
        local_elemental_stiffness(1, 5) = tmp6z;
        local_elemental_stiffness(1, 7) = -tmp12z;
        local_elemental_stiffness(1, 11) = tmp6z;

        local_elemental_stiffness(2, 2) = tmp12y;
        local_elemental_stiffness(2, 4) = -tmp6y;
        local_elemental_stiffness(2, 8) = -tmp12y;
        local_elemental_stiffness(2, 10) = -tmp6y;

        local_elemental_stiffness(3, 3) = tmpGJ;
        local_elemental_stiffness(3, 9) = -tmpGJ;

        local_elemental_stiffness(4, 2) = -tmp6y;
        local_elemental_stiffness(4, 4) = 4.0 * tmp1y;
        local_elemental_stiffness(4, 8) = tmp6y;
        local_elemental_stiffness(4, 10) = 2.0 * tmp1y;

        local_elemental_stiffness(5, 1) = tmp6z;
        local_elemental_stiffness(5, 5) = 4.0 * tmp1z;
        local_elemental_stiffness(5, 7) = -tmp6z;
        local_elemental_stiffness(5, 11) = 2.0 * tmp1z;

        local_elemental_stiffness(6, 0) = -tmpEA;
        local_elemental_stiffness(6, 6) = tmpEA;

        local_elemental_stiffness(7, 1) = -tmp12z;
        local_elemental_stiffness(7, 5) = -tmp6z;
        local_elemental_stiffness(7, 7) = tmp12z;
        local_elemental_stiffness(7, 11) = -tmp6z;

        local_elemental_stiffness(8, 2) = -tmp12y;
        local_elemental_stiffness(8, 4) = tmp6y;
        local_elemental_stiffness(8, 8) = tmp12y;
        local_elemental_stiffness(8, 10) = tmp6y;

        local_elemental_stiffness(9, 3) = -tmpGJ;
        local_elemental_stiffness(9, 9) = tmpGJ;

        local_elemental_stiffness(10, 2) = -tmp6y;
        local_elemental_stiffness(10, 4) = 2.0 * tmp1y;
        local_elemental_stiffness(10, 8) = tmp6y;
        local_elemental_stiffness(10, 10) = 4.0 * tmp1y;

        local_elemental_stiffness(11, 1) = tmp6z;
        local_elemental_stiffness(11, 5) = 2.0 * tmp1z;
        local_elemental_stiffness(11, 7) = -tmp6z;
        local_elemental_stiffness(11, 11) = 4.0 * tmp1z;

        update_rotation(nodes, *this, rotation, rotation_transposed);
        return rotation_transposed * local_elemental_stiffness * rotation;
    }

    LocalMatrix EulerBernoulliBeamElement::calculate_inv_mass_matrix(const std::vector<Node> &nodes) {
        // store node indices of current element
        const int nn1 = node_numbers[0];
        const int nn2 = node_numbers[1];

        // calculate the length of the element
        const Node dn = nodes[nn1] - nodes[nn2];
        const double length = dn.norm();

        const double mass = props.density * length * props.area;

        // store the entries in the (local) elemental stiffness matrix as temporary values to avoid recalculation
        const double tmp2 = 2.0 / mass;
        const double tmp4 = 4.0 / mass;
        const double tmp16 = 16.0 / mass;
        const double tmp60 = 60 / (mass * length);
        const double tmp120 = 120.0 / (mass * length);
        const double tmp840 = 840 / (mass * length * length);
        const double tmp1200 = 1200.0 / (mass * length * length);

        LocalMatrix inverse_local_elemental_mass, rotation, rotation_transposed;
        inverse_local_elemental_mass.setZero();

        // update local elemental mass matrix
        inverse_local_elemental_mass(0, 0) = tmp4;
        inverse_local_elemental_mass(0, 6) = -tmp2;

        inverse_local_elemental_mass(1, 1) = tmp16;
        inverse_local_elemental_mass(1, 5) = -tmp120;
        inverse_local_elemental_mass(1, 7) = -tmp4;
        inverse_local_elemental_mass(1, 11) = -tmp60;

        inverse_local_elemental_mass(2, 2) = tmp16;
        inverse_local_elemental_mass(2, 4) = -tmp120;
        inverse_local_elemental_mass(2, 8) = -tmp4;
        inverse_local_elemental_mass(2, 10) = -tmp60;

        inverse_local_elemental_mass(3, 3) = tmp4;
        inverse_local_elemental_mass(3, 9) = -tmp2;

        inverse_local_elemental_mass(4, 2) = -tmp120;
        inverse_local_elemental_mass(4, 4) = tmp1200;
        inverse_local_elemental_mass(4, 8) = tmp60;
        inverse_local_elemental_mass(4, 10) = tmp840;

        inverse_local_elemental_mass(5, 1) = -tmp120;
        inverse_local_elemental_mass(5, 5) = tmp1200;
        inverse_local_elemental_mass(5, 7) = tmp60;
        inverse_local_elemental_mass(5, 11) = tmp840;

        inverse_local_elemental_mass(6, 0) = -tmp2;
        inverse_local_elemental_mass(6, 6) = tmp4;

        inverse_local_elemental_mass(7, 1) = -tmp4;
        inverse_local_elemental_mass(7, 5) = tmp60;
        inverse_local_elemental_mass(7, 7) = tmp16;
        inverse_local_elemental_mass(7, 11) = tmp120;

        inverse_local_elemental_mass(8, 2) = -tmp4;
        inverse_local_elemental_mass(8, 4) = tmp60;
        inverse_local_elemental_mass(8, 8) = tmp16;
        inverse_local_elemental_mass(8, 10) = tmp120;

        inverse_local_elemental_mass(9, 3) = -tmp2;
        inverse_local_elemental_mass(9, 9) = tmp4;

        inverse_local_elemental_mass(10, 2) = -tmp60;
        inverse_local_elemental_mass(10, 4) = tmp840;
        inverse_local_elemental_mass(10, 8) = tmp120;
        inverse_local_elemental_mass(10, 10) = tmp1200;

        inverse_local_elemental_mass(11, 1) = -tmp60;
        inverse_local_elemental_mass(11, 5) = tmp840;
        inverse_local_elemental_mass(11, 7) = tmp120;
        inverse_local_elemental_mass(11, 11) = tmp1200;

        update_rotation(nodes, *this, rotation, rotation_transposed);
        return rotation_transposed * inverse_local_elemental_mass * rotation;
    }

} // namespace explicit_fea