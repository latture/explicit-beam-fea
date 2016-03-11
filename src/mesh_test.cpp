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

#include <gtest/gtest.h>
#include "mesh.h"
#include "euler_bernoulli_beam_element.h"

namespace explicit_fea {
    TEST(Mesh, AssemblesGlobalStiffness) {

        double youngs_modulus = 10.0;
        double shear_modulus = 10.0;
        double area = 1.0;
        double Iy = 1.0;
        double Iz = 1.0;
        double J = 1.0;
        double density = 1.0;
        Eigen::Vector3d ny;
        ny << 0.0, 1.0, 0.0;

        Props props1(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, ny);
        Props props2(youngs_modulus/10., shear_modulus/10., area, Iz, Iy, J, density, ny);

        std::vector<Node> nodes = {Node(0.0, 0.0, 0.0), Node(1.0, 0.0, 0.0), Node(2.0, 0.0, 0.0), Node(2.0, 0.0, 1.0)};

        std::unique_ptr<EulerBernoulliBeamElement> init[] = {std::unique_ptr<EulerBernoulliBeamElement>(new EulerBernoulliBeamElement(0, 1, props1)),
                                                             std::unique_ptr<EulerBernoulliBeamElement>(new EulerBernoulliBeamElement(1, 2, props1)),
                                                             std::unique_ptr<EulerBernoulliBeamElement>(new EulerBernoulliBeamElement(2, 3, props2))};

        std::vector<std::unique_ptr<BeamElement>> elems = {std::make_move_iterator(std::begin(init)), std::make_move_iterator(std::end(init))};
        BCList bcs;

        Mesh mesh(nodes, elems, std::move(bcs));

        SparseMatrix actual_sp = mesh.get_global_stiffness_matrix();
        Eigen::MatrixXd actual = Eigen::MatrixXd(actual_sp);

        std::vector<std::vector<double> > expected =
               {{10.,0.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,120.,0.,0.,0.,60.,0.,-120.,0.,0.,0.,60.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,0.,120.,0.,-60.,0.,0.,0.,-120.,0.,-60.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,0.,0.,10.,0.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,0.,-60.,0.,40.,0.,0.,0.,60.,0.,20.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,60.,0.,0.,0.,40.,0.,-60.,0.,0.,0.,20.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {-10.,0.,0.,0.,0.,0.,20.,0.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,-120.,0.,0.,0.,-60.,0.,240.,0.,0.,0.,0.,0.,-120.,0.,0.,0.,60.,0.,0.,0.,0.,0.,0.},
                {0.,0.,-120.,0.,60.,0.,0.,0.,240.,0.,0.,0.,0.,0.,-120.,0.,-60.,0.,0.,0.,0.,0.,0.,0.},
                {0.,0.,0.,-10.,0.,0.,0.,0.,0.,20.,0.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,0.,0.,0.},
                {0.,0.,-60.,0.,20.,0.,0.,0.,0.,0.,80.,0.,0.,0.,60.,0.,20.,0.,0.,0.,0.,0.,0.,0.},
                {0.,60.,0.,0.,0.,20.,0.,0.,0.,0.,0.,80.,0.,-60.,0.,0.,0.,20.,0.,0.,0.,0.,0.,0.},
                {0.,0.,0.,0.,0.,0.,-10.,0.,0.,0.,0.,0.,22.,0.,0.,0.,6.,0.,-12.,0.,0.,0.,6.,0.},
                {0.,0.,0.,0.,0.,0.,0.,-120.,0.,0.,0.,-60.,0.,132.,0.,-6.,0.,-60.,0.,-12.,0.,-6.,0.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,-120.,0.,60.,0.,0.,0.,121.,0.,60.,0.,0.,0.,-1.,0.,0.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,-10.,0.,0.,0.,-6.,0.,14.,0.,0.,0.,6.,0.,2.,0.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,-60.,0.,20.,0.,6.,0.,60.,0.,44.,0.,-6.,0.,0.,0.,2.,0.},
                {0.,0.,0.,0.,0.,0.,0.,60.,0.,0.,0.,20.,0.,-60.,0.,0.,0.,41.,0.,0.,0.,0.,0.,-1.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-12.,0.,0.,0.,-6.,0.,12.,0.,0.,0.,-6.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-12.,0.,6.,0.,0.,0.,12.,0.,6.,0.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,0.,0.,0.,0.,0.,1.,0.,0.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-6.,0.,2.,0.,0.,0.,6.,0.,4.,0.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,6.,0.,0.,0.,2.,0.,-6.,0.,0.,0.,4.,0.},
                {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-1.,0.,0.,0.,0.,0.,1.}};

       for (size_t i = 0; i < expected.size(); i++) {
           for (size_t j = 0; j < expected[i].size(); j++) {
               EXPECT_DOUBLE_EQ(expected[i][j], actual(i, j)) << "\ti = " << i << "\tj =  " << j;
           }
       }
    }

    TEST(Mesh, AppliesBC) {
        double youngs_modulus = 10.0;
        double shear_modulus = 10.0;
        double area = 1.0;
        double Iy = 1.0;
        double Iz = 1.0;
        double J = 1.0;
        double density = 1.0;
        Eigen::Vector3d ny;
        ny << 0.0, 1.0, 0.0;

        Props props1(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, ny);
        Props props2(youngs_modulus/10., shear_modulus/10., area, Iz, Iy, J, density, ny);

        std::vector<Node> nodes = {Node(0.0, 0.0, 0.0), Node(1.0, 0.0, 0.0), Node(2.0, 0.0, 0.0), Node(2.0, 0.0, 1.0)};

        std::unique_ptr<EulerBernoulliBeamElement> init[] = {std::unique_ptr<EulerBernoulliBeamElement>(new EulerBernoulliBeamElement(0, 1, props1)),
                                                             std::unique_ptr<EulerBernoulliBeamElement>(new EulerBernoulliBeamElement(1, 2, props1)),
                                                             std::unique_ptr<EulerBernoulliBeamElement>(new EulerBernoulliBeamElement(2, 3, props2))};

        std::vector<std::unique_ptr<BeamElement>> elems = {std::make_move_iterator(std::begin(init)), std::make_move_iterator(std::end(init))};
        BCList bcs;
        bcs.push_back(std::unique_ptr<ConstantBC>(new ConstantBC(0, 0, 0.0, BC::Type::DISPLACEMENT)));
        bcs.push_back(std::unique_ptr<ConstantBC>(new ConstantBC(1, 4, 0.0, BC::Type::DISPLACEMENT)));

        Mesh mesh(nodes, elems, std::move(bcs));

        SparseMatrix actual_sp = mesh.get_inv_mass_matrix();
        Eigen::MatrixXd actual = Eigen::MatrixXd(actual_sp);

        EXPECT_EQ(actual.rows(), actual.cols()) << "Inverse mass matrix must be square.";

        int global_idx;
        for (BCList::const_iterator it = mesh.get_bcs().begin(); it != mesh.get_bcs().end(); ++it) {
            global_idx = DOF::NUM_DOFS * (*it)->node + (*it)->dof;
            for (size_t i = 0; i < actual.rows(); ++i) {
                if (i == global_idx) {
                    EXPECT_DOUBLE_EQ(1.0, actual(global_idx, global_idx));
                }
                else {
                    EXPECT_DOUBLE_EQ(0.0, actual(i, global_idx));
                    EXPECT_DOUBLE_EQ(0.0, actual(global_idx, i));
                }
            }
        }
    }
} // namespace explicit_fea



