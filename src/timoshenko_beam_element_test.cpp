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
#include "timoshenko_beam_element.h"
#include "fe_utils.h"
#include <memory>

using namespace explicit_fea;

TEST(TimoshenkoBeamElement, AssemblesElementalStiffness) {
    Node n1, n2;
    n1 << 0.0, 0.0, 0.0;
    n2 << 1.0, 0.0, 0.0;
    std::vector<Node> nodes = {n1, n2};

    double youngs_modulus = 10.0;
    double shear_modulus = 10.0;
    double area = 1.0;
    double Iy = 1.0;
    double Iz = 1.0;
    double J = 1.0;
    double density = 1.0;
    Eigen::Vector3d ny;
    ny << 0.0, 0.0, 1.0;
    Props p(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, ny);

    TimoshenkoBeamElement elem(0, 1, p);

    LocalMatrix actual = elem.calculate_stiffness_matrix(nodes);

    LocalMatrix expected, rotation, rotation_transposed;
    update_rotation(nodes, elem, rotation, rotation_transposed);
    expected << 10, 0, 0, 0, 0, 0, -10, 0, 0, 0, 0, 0,
                 0, 120./13., 0, 0, 0, 60./13., 0, -120./13., 0, 0, 0, 60./13.,
                 0, 0, 120./13., 0, 60./13., 0, 0, 0, -120./13., 0, 60./13., 0,
                 0, 0, 0, 10, 0, 0, 0, 0, 0, -10, 0, 0,
                 0, 0, 60./13., 0, 160./13., 0, 0, 0, -60./13., 0, -100./13., 0,
                 0, 60./13., 0, 0, 0, 160./13., 0, -60./13., 0, 0, 0, -100./13.,
                 -10., 0, 0, 0, 0, 0, 10., 0, 0, 0, 0, 0,
                 0, -120./13., 0, 0, 0, -60./13., 0, 120./13., 0, 0, 0, -60./13.,
                 0, 0, -120./13., 0, -60./13., 0, 0, 0, 120./13., 0, -60./13., 0,
                 0, 0, 0, -10, 0, 0, 0, 0, 0, 10, 0, 0,
                 0, 0, 60./13., 0, -100./13., 0, 0, 0, -60./13., 0, 160./13., 0,
                 0, 60./13., 0, 0, 0, -100./13., 0, -60./13., 0, 0, 0, 160./13.;

    // transform expected output from local to global coordinates
    expected = rotation_transposed * expected * rotation;

    for (size_t i = 0; i < expected.rows(); i++) {
        for (size_t j = 0; j < expected.cols(); j++) {
            EXPECT_DOUBLE_EQ(expected(i, j), actual(i, j)) << "\ti = " << i << "\tj =  " << j;
        }
    }
}

TEST(TimoshenkoBeamElement, AssemblesTranslatedElementalStiffness) {
    Node n1, n2;
    n1 << 1.0, 0.0, 0.0;
    n2 << 2.0, 0.0, 0.0;
    std::vector<Node> nodes = {n1, n2};

    double youngs_modulus = 10.0;
    double shear_modulus = 10.0;
    double area = 1.0;
    double Iy = 1.0;
    double Iz = 1.0;
    double J = 1.0;
    double density = 1.0;
    Eigen::Vector3d ny;
    ny << 0.0, 0.0, 1.0;
    Props p(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, ny);

    TimoshenkoBeamElement elem(0, 1, p);

    LocalMatrix actual = elem.calculate_stiffness_matrix(nodes);

    LocalMatrix expected, rotation, rotation_transposed;
    update_rotation(nodes, elem, rotation, rotation_transposed);
    expected << 10, 0, 0, 0, 0, 0, -10, 0, 0, 0, 0, 0,
            0, 120./13., 0, 0, 0, 60./13., 0, -120./13., 0, 0, 0, 60./13.,
            0, 0, 120./13., 0, 60./13., 0, 0, 0, -120./13., 0, 60./13., 0,
            0, 0, 0, 10, 0, 0, 0, 0, 0, -10, 0, 0,
            0, 0, 60./13., 0, 160./13., 0, 0, 0, -60./13., 0, -100./13., 0,
            0, 60./13., 0, 0, 0, 160./13., 0, -60./13., 0, 0, 0, -100./13.,
            -10., 0, 0, 0, 0, 0, 10., 0, 0, 0, 0, 0,
            0, -120./13., 0, 0, 0, -60./13., 0, 120./13., 0, 0, 0, -60./13.,
            0, 0, -120./13., 0, -60./13., 0, 0, 0, 120./13., 0, -60./13., 0,
            0, 0, 0, -10, 0, 0, 0, 0, 0, 10, 0, 0,
            0, 0, 60./13., 0, -100./13., 0, 0, 0, -60./13., 0, 160./13., 0,
            0, 60./13., 0, 0, 0, -100./13., 0, -60./13., 0, 0, 0, 160./13.;

    // transform expected output from local to global coordinates
    expected = rotation_transposed * expected * rotation;

    for (size_t i = 0; i < expected.rows(); i++) {
        for (size_t j = 0; j < expected.cols(); j++) {
            EXPECT_DOUBLE_EQ(expected(i, j), actual(i, j)) << "\ti = " << i << "\tj =  " << j;
        }
    }
}

TEST(TimoshenkoBeamElement, AssemblesElementalInvMass) {
    Node n1, n2;
    n1 << 0.0, 0.0, 0.0;
    n2 << 1.0, 0.0, 0.0;
    std::vector<Node> nodes = {n1, n2};

    double youngs_modulus = 10.0;
    double shear_modulus = 10.0;
    double area = 1.0;
    double Iy = 1.0;
    double Iz = 1.0;
    double J = 1.0;
    double density = 1.0;
    Eigen::Vector3d ny;
    ny << 0.0, 0.0, 1.0;

    Props p(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, ny);
    TimoshenkoBeamElement elem(0, 1, p);

    LocalMatrix actual = elem.calculate_inv_mass_matrix(nodes);

    LocalMatrix expected, rotation, rotation_transposed;
    update_rotation(nodes, elem, rotation, rotation_transposed);
    expected << 4, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0,
            0, 2704./11809., 0, 0, 0, -225600./11809., 0, -2116./11809., 0, 0, 0, -222660./11809.,
            0, 0, 2704./11809., 0, -225600./11809., 0, 0, 0, -2116./11809., 0, -222660./11809., 0,
            0, 0, 0, 4, 0, 0, 0, 0, 0, -2, 0, 0,
            0, 0, -225600./11809., 0, 3943349040./1995721., 0, 0, 0, 222660./11809., 0, 3940156200./1995721., 0,
            0, -225600./11809., 0, 0, 0, 3943349040./1995721., 0, 222660./11809., 0, 0, 0, 3940156200./1995721.,
            -2, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0,
            0, -2116./11809., 0, 0, 0, 222660./11809., 0, 2704./11809., 0, 0, 0, 225600./11809.,
            0, 0, -2116./11809., 0, 222660./11809., 0, 0, 0, 2704./11809., 0, 225600./11809., 0,
            0, 0, 0, -2, 0, 0, 0, 0, 0, 4, 0, 0,
            0, 0, -222660./11809., 0, 3940156200./1995721., 0, 0, 0, 225600./11809., 0, 3943349040./1995721., 0,
            0, -222660./11809., 0, 0, 0, 3940156200./1995721., 0, 225600./11809., 0, 0, 0, 3943349040./1995721.;

    // transform expected output from local to global coordinates
    expected = rotation_transposed * expected * rotation;

    for (size_t i = 0; i < expected.rows(); i++) {
        for (size_t j = 0; j < expected.cols(); j++) {
            EXPECT_DOUBLE_EQ(expected(i, j), actual(i, j)) << "\ti = " << i << "\tj =  " << j;
        }
    }
}
