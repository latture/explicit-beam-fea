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
#include "explicit_system.h"

using namespace explicit_fea;

TEST(ExplicitSystem, IntegratesSingleElem)
{
    double youngs_modulus = 200.0e9;
    double shear_modulus = 80.0e9;
    double area = 0.0314159265358979;
    double Iy = 0.0000785398;
    double Iz = 0.0000785398;
    double J = 2.0 * Iy;
    double density = 7800.0;
    Eigen::Vector3d ny;
    ny << 0.0, 1.0, 0.0;
    Props props(youngs_modulus, shear_modulus, area, Iz, Iy, J, density, ny);
    std::vector<Node> nodes = {Node(0.0, 0.0, 0.0), Node(1.0, 0.0, 0.0)};

    std::unique_ptr<TimoshenkoBeamElement> init[] = {std::unique_ptr<TimoshenkoBeamElement>(new TimoshenkoBeamElement(0, 1, props))};
    std::vector<std::unique_ptr<BeamElement>> elems = {std::make_move_iterator(std::begin(init)), std::make_move_iterator(std::end(init))};

    BCList bcs;
    for (unsigned int i = 0; i < DOF::NUM_DOFS; ++i) {
        bcs.push_back(std::unique_ptr<ConstantBC>(new ConstantBC(0, i, 0.0, BC::Type::DISPLACEMENT)));
    }
    double velocity = 0.001;
    bcs.push_back(std::unique_ptr<ConstantBC>(new ConstantBC(1, DOF::DISPLACEMENT_X, velocity, BC::Type::VELOCITY)));

    Mesh mesh(nodes, elems, std::move(bcs));
    ForceList forces;

    ColumnVector initial_displacements(mesh.get_global_stiffness_matrix().cols());
    ColumnVector initial_velocities(mesh.get_global_stiffness_matrix().cols());
    initial_displacements.setZero();
    initial_velocities.setZero();

    ExplicitSystem::Options options;
    double t_0 = 0.0;
    std::unique_ptr<ExplicitSystem> es (new ExplicitSystem(std::move(mesh),
                                                           std::move(forces),
                                                           initial_displacements,
                                                           initial_velocities,
                                                           t_0,
                                                           options));

    double timestep = estimate_stable_timestep(nodes, elems);
    int num_iter = 1000;
    for (int j = 0; j < num_iter; ++j) {
        es->update(timestep);
    }

	// Check displacements
    ColumnVector actual_displacements = es->getDisplacements();
    ColumnVector expected_displacements(ColumnVector::Zero(mesh.get_global_stiffness_matrix().cols()));
    expected_displacements[6] = num_iter * timestep * velocity;
    for (int i = 0; i < expected_displacements.size(); ++i) 
        EXPECT_FLOAT_EQ(expected_displacements[i], actual_displacements[i]);
    

	// Check forces
	ColumnVector actual_forces = es->getForces();
	ColumnVector expected_forces = (ColumnVector::Zero(mesh.get_global_stiffness_matrix().cols()));
	double strain = num_iter * timestep * velocity;
	double stress = strain * youngs_modulus;
	double reaction_force = stress * area;
	expected_forces[0] = -reaction_force;
	expected_forces[6] = reaction_force;

	for (size_t i = 0; i < expected_forces.size(); i++)
		EXPECT_FLOAT_EQ(expected_forces[i], actual_forces[i]);
}


