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

#include "timoshenko_beam_element.h"
#include "fe_utils.h"

namespace explicit_fea {


    TimoshenkoBeamElement::TimoshenkoBeamElement(int nn1, int nn2, const Props& props) : BeamElement(nn1, nn2, props) {
    }

    LocalMatrix TimoshenkoBeamElement::calculate_stiffness_matrix(const std::vector <Node> &nodes) {
        // store node indices of current element
        const int nn1 = node_numbers[0];
        const int nn2 = node_numbers[1];

        const Node dn = nodes[nn1] - nodes[nn2];
        const double length = dn.norm();
        const double EA = props.youngs_modulus * props.area;
        const double EIz = props.youngs_modulus * props.Iz;
        const double EIy = props.youngs_modulus * props.Iy;
        const double GJ = props.shear_modulus * props.J;
        const double phiy = 12.0 * EIy / (props.shear_modulus * props.area * length * length);
        const double phiz = 12.0 * EIz / (props.shear_modulus * props.area * length * length);

        // store the entries in the (local) elemental stiffness matrix as temporary values to avoid recalculation
        const double tmpEA = EA / length;
        const double tmpGJ = GJ / length;

        const double tmp12z = 12.0 * EIz / ((length * length * length) * (1.0 + phiz));
        const double tmp6z = 6.0 * EIz / ((length * length) * (1.0 + phiz));
        const double tmp1z2 = EIz * (2.0 - phiz) / (length * (1.0 + phiz));
        const double tmp1z4 = EIz * (4.0 + phiz) / (length * (1.0 + phiz));

        const double tmp12y = 12.0 * EIy / ((length * length * length) * (1.0 + phiy));
        const double tmp6y = 6.0 * EIy / ((length * length) * (1.0 + phiy));
        const double tmp1y2 = EIy * (2.0 - phiy) / (length * (1.0 + phiy));
        const double tmp1y4 = EIy * (4.0 + phiy) / (length * (1.0 + phiy));

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
        local_elemental_stiffness(4, 4) = tmp1y4;
        local_elemental_stiffness(4, 8) = tmp6y;
        local_elemental_stiffness(4, 10) = tmp1y2;

        local_elemental_stiffness(5, 1) = tmp6z;
        local_elemental_stiffness(5, 5) = tmp1z4;
        local_elemental_stiffness(5, 7) = -tmp6z;
        local_elemental_stiffness(5, 11) = tmp1z2;

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
        local_elemental_stiffness(10, 4) = tmp1y2;
        local_elemental_stiffness(10, 8) = tmp6y;
        local_elemental_stiffness(10, 10) = tmp1y4;

        local_elemental_stiffness(11, 1) = tmp6z;
        local_elemental_stiffness(11, 5) = tmp1z2;
        local_elemental_stiffness(11, 7) = -tmp6z;
        local_elemental_stiffness(11, 11) = tmp1z4;

        update_rotation(nodes, *this, rotation, rotation_transposed);
        return rotation_transposed * local_elemental_stiffness * rotation;
    }

    LocalMatrix TimoshenkoBeamElement::calculate_inv_mass_matrix(const std::vector<Node> &nodes) {
        // store node indices of current element
        const int nn1 = node_numbers[0];
        const int nn2 = node_numbers[1];

        const Node dn = nodes[nn1] - nodes[nn2];
        const double length = dn.norm();
        const double mass = props.density * length * props.area;
        const double EIz = props.youngs_modulus * props.Iz;
        const double EIy = props.youngs_modulus * props.Iy;
        const double phiy = 12.0 * EIy / (props.shear_modulus * props.area * length * length);
        const double phiz = 12.0 * EIz / (props.shear_modulus * props.area * length * length);


        // store the entries in the (local) elemental stiffness matrix as temporary values to avoid recalculation
        const double tmp2 = 2.0 / mass;
        const double tmp4 = 4.0 / mass;

        const double tmpd1z = mass * (6.0 + phiz * (12.0 + phiz)) * (2.0 + phiz * (4.0 + 3.0 * phiz));
        const double invtmpd2z = 1.0 / (pow(length, 2) * mass * pow(1.0 + phiz, 2) * (6.0 + phiz * (12.0 + phiz)) * (2.0 + phiz * (4.0 + 3.0 * phiz)));
        const double tmp192z = (192.0 * pow(1.0 + phiz, 2)) / tmpd1z;
        const double tmp24z = (24.0 * (2.0 + phiz * (4.0 + 7.0 * phiz))) / tmpd1z;
        const double tmp6024z = (60.0 * (24.0 + phiz * (62.0 + 7.0 * phiz * (8.0 + 3.0 * phiz)))) / (length * tmpd1z);
        const double tmp6012z = (60.0 * (12.0 + phiz * (38.0 + 3.0 * phiz * (18.0 + 7.0 * phiz)))) / (length * tmpd1z);
        const double tmp30480z = (30.0 * (480.0 + phiz * (2592.0 + phiz * (5928.0 + phiz * (7428.0 + phiz * (5350.0 + 21.0 * phiz * (98.0 + 15.0 * phiz))))))) * invtmpd2z;
        const double tmp30336z = (30.0 * (336.0 + phiz * (2016.0 + phiz * (5172.0 + phiz * (7068.0 + phiz * (5324.0 + 21.0 * phiz * (98.0 + 15.0 * phiz))))))) * invtmpd2z;

        const double tmpd1y = mass * (6.0 + phiy * (12.0 + phiy)) * (2.0 + phiy * (4.0 + 3.0 * phiy));
        const double invtmpd2y = 1.0 / (pow(length, 2) * mass * pow(1.0 + phiy, 2) * (6.0 + phiy * (12.0 + phiy)) * (2.0 + phiy * (4.0 + 3.0 * phiy)));
        const double tmp192y = (192.0 * pow(1.0 + phiy, 2)) / tmpd1y;
        const double tmp24y = (24.0 * (2.0 + phiy * (4.0 + 7.0 * phiy))) / tmpd1y;
        const double tmp6024y = (60.0 * (24.0 + phiy * (62.0 + 7.0 * phiy * (8.0 + 3.0 * phiy)))) / (length * tmpd1y);
        const double tmp6012y = (60.0 * (12.0 + phiy * (38.0 + 3.0 * phiy * (18.0 + 7.0 * phiy)))) / (length * tmpd1y);
        const double tmp30480y = (30.0 * (480.0 + phiy * (2592.0 + phiy * (5928.0 + phiy * (7428.0 + phiy * (5350.0 + 21.0 * phiy * (98.0 + 15.0 * phiy))))))) * invtmpd2y;
        const double tmp30336y = (30.0 * (336.0 + phiy * (2016.0 + phiy * (5172.0 + phiy * (7068.0 + phiy * (5324.0 + 21.0 * phiy * (98.0 + 15.0 * phiy))))))) * invtmpd2y;


        LocalMatrix inverse_local_elemental_mass, rotation, rotation_transposed;
        inverse_local_elemental_mass.setZero();

        // update local elemental mass matrix
        inverse_local_elemental_mass(0, 0) = tmp4;
        inverse_local_elemental_mass(0, 6) = -tmp2;

        inverse_local_elemental_mass(1, 1) = tmp192z;
        inverse_local_elemental_mass(1, 5) = -tmp6024z;
        inverse_local_elemental_mass(1, 7) = -tmp24z;
        inverse_local_elemental_mass(1, 11) = -tmp6012z;

        inverse_local_elemental_mass(2, 2) = tmp192y;
        inverse_local_elemental_mass(2, 4) = -tmp6024y;
        inverse_local_elemental_mass(2, 8) = -tmp24y;
        inverse_local_elemental_mass(2, 10) = -tmp6012y;

        inverse_local_elemental_mass(3, 3) = tmp4;
        inverse_local_elemental_mass(3, 9) = -tmp2;

        inverse_local_elemental_mass(4, 2) = -tmp6024y;
        inverse_local_elemental_mass(4, 4) = tmp30480y;
        inverse_local_elemental_mass(4, 8) = tmp6012y;
        inverse_local_elemental_mass(4, 10) = tmp30336y;

        inverse_local_elemental_mass(5, 1) = -tmp6024z;
        inverse_local_elemental_mass(5, 5) = tmp30480z;
        inverse_local_elemental_mass(5, 7) = tmp6012z;
        inverse_local_elemental_mass(5, 11) = tmp30336z;

        inverse_local_elemental_mass(6, 0) = -tmp2;
        inverse_local_elemental_mass(6, 6) = tmp4;

        inverse_local_elemental_mass(7, 1) = -tmp24z;
        inverse_local_elemental_mass(7, 5) = tmp6012z;
        inverse_local_elemental_mass(7, 7) = tmp192z;
        inverse_local_elemental_mass(7, 11) = tmp6024z;

        inverse_local_elemental_mass(8, 2) = -tmp24y;
        inverse_local_elemental_mass(8, 4) = tmp6012y;
        inverse_local_elemental_mass(8, 8) = tmp192y;
        inverse_local_elemental_mass(8, 10) = tmp6024y;

        inverse_local_elemental_mass(9, 3) = -tmp2;
        inverse_local_elemental_mass(9, 9) = tmp4;

        inverse_local_elemental_mass(10, 2) = -tmp6012y;
        inverse_local_elemental_mass(10, 4) = tmp30336y;
        inverse_local_elemental_mass(10, 8) = tmp6024y;
        inverse_local_elemental_mass(10, 10) = tmp30480y;

        inverse_local_elemental_mass(11, 1) = -tmp6012z;
        inverse_local_elemental_mass(11, 5) = tmp30336z;
        inverse_local_elemental_mass(11, 7) = tmp6024z;
        inverse_local_elemental_mass(11, 11) = tmp30480z;

        update_rotation(nodes, *this, rotation, rotation_transposed);
        return rotation_transposed * inverse_local_elemental_mass * rotation;
    }


} // namespace explicit_fea