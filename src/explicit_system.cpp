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

#include "explicit_system.h"
#include <boost/format.hpp>
#include <fstream>

namespace explicit_fea {

    ExplicitSystem::ExplicitSystem(Mesh _mesh,
                                   ForceList _external_forces,
                                   const ColumnVector& _initial_displacements,
                                   const ColumnVector& _initial_velocities,
                                   double _t_0,
                                   const Options& _options) :
        mesh(std::move(_mesh)),
        external_forces(std::move(_external_forces)),
        displacements_0(_initial_displacements),
        velocities_0(_initial_velocities),
        velocities_1(ColumnVector::Zero(_initial_displacements.size())),
        accelerations_0(ColumnVector::Zero(_initial_displacements.size())),
        accelerations_1(ColumnVector::Zero(_initial_displacements.size())),
        forces(ColumnVector::Zero(_initial_displacements.size())),
        RHS(ColumnVector::Zero(_initial_displacements.size())),
        t_0(_t_0),
        dt_0(0.0),
        options(_options) {
        if (displacements_0.size() != velocities_0.size())
            throw std::runtime_error("Size of initial velocities and initial displacements are not equal.");
        if (displacements_0.size() != mesh.get_global_stiffness_matrix().cols())
            throw std::runtime_error("Size of displacements and velocities does not match the number of columns in the global matrices.");

        LHS.resize(mesh.get_global_stiffness_matrix().rows(), mesh.get_global_stiffness_matrix().cols());
        assemble_damping_matrix();
        apply_bcs(t_0);
    }

    void ExplicitSystem::Options::load(const rapidjson::Document &config_doc) {
        if (config_doc.HasMember("options")) {
            if (config_doc["options"].HasMember("beta")) {
                if (!config_doc["options"]["beta"].IsNumber()) {
                    throw std::runtime_error("beta provided in options configuration is not a number.");
                }
                beta = config_doc["options"]["beta"].GetDouble();
            }
            if (config_doc["options"].HasMember("gamma")) {
                if (!config_doc["options"]["gamma"].IsNumber()) {
                    throw std::runtime_error("gamma provided in options configuration is not a number.");
                }
                gamma = config_doc["options"]["gamma"].GetDouble();
            }
            if (config_doc["options"].HasMember("damping_beta")) {
                if (!config_doc["options"]["damping_beta"].IsNumber()) {
                    throw std::runtime_error("damping_beta provided in options configuration is not a number.");
                }
                damping_beta = config_doc["options"]["damping_beta"].GetDouble();
            }
            if (config_doc["options"].HasMember("damping_alpha")) {
                if (!config_doc["options"]["damping_alpha"].IsNumber()) {
                    throw std::runtime_error("damping_alpha provided in options configuration is not a number.");
                }
                damping_alpha = config_doc["options"]["damping_alpha"].GetDouble();
            }
        }
    }

    void ExplicitSystem::update(double dt_1) {
        double t_1 = t_0 + dt_1;
        if (!compare.equal(dt_1, dt_0)) {
            assemble_LHS(dt_1);
        }

        apply_external_forces(t_1);
        RHS = forces;
        RHS.noalias() -= damping_matrix * (velocities_0 + (1.0 - options.gamma) * dt_1 * accelerations_0);
        RHS.noalias() -= mesh.get_global_stiffness_matrix() * (displacements_0 + dt_1 * velocities_0 + (0.5 - options.beta) * dt_1 * dt_1 * accelerations_0);
        apply_bcs(t_1);

        accelerations_1 = solver.solve(RHS);
        velocities_1.noalias() += (1.0 - options.gamma) * dt_1 * accelerations_0 + options.gamma * dt_1 * accelerations_1;

        // advance all current time values to next iteration
        displacements_0.noalias() += dt_1 * velocities_0 + dt_1 * dt_1 * (0.5 - options.beta) * accelerations_0 + dt_1 * dt_1 * options.beta * accelerations_1;
        velocities_0 = velocities_1;
        accelerations_0 = accelerations_1;
        t_0 = t_1;
        dt_0 = dt_1;
    }

	Mesh const& ExplicitSystem::getMesh() const {
		return mesh;
	}

    ColumnVector const& ExplicitSystem::getDisplacements() const {
        return displacements_0;
    }

    ColumnVector ExplicitSystem::getForces() const {
        ColumnVector out(displacements_0.size());
        out = mesh.get_global_stiffness_matrix() * displacements_0;
        out += mesh.get_mass_matrix() * accelerations_0;
        return out;
    }

    ColumnVector const& ExplicitSystem::getVelocities() const {
        return velocities_0;
    }

    double ExplicitSystem::getTime() const {
        return t_0;
    }

    void ExplicitSystem::assemble_damping_matrix() {
        damping_matrix.resize(mesh.get_global_stiffness_matrix().rows(), mesh.get_global_stiffness_matrix().cols());
        damping_matrix = options.damping_alpha * mesh.get_mass_matrix() + options.damping_beta * mesh.get_global_stiffness_matrix();
        damping_matrix.prune(1.e-14);
        damping_matrix.makeCompressed();
    }

    void ExplicitSystem::assemble_LHS(double dt){
        LHS = mesh.get_mass_matrix() + options.gamma * dt * damping_matrix + options.beta * dt * dt * mesh.get_global_stiffness_matrix();
        LHS.prune(1.e-14);
        LHS.makeCompressed();
        solver.compute(LHS);
    }

    void ExplicitSystem::apply_bcs(double time) {
        for (BCList::const_iterator it = mesh.get_bcs().begin(); it != mesh.get_bcs().end(); ++it) {
            RHS[(*it)->global_index] = 0.0;
            if ((*it)->type == BC::Type::DISPLACEMENT) {
                displacements_0[(*it)->global_index] = (*it)->get_value(time);
                velocities_0[(*it)->global_index] = 0.0;
            } else if ((*it)->type == BC::Type::VELOCITY) {
                velocities_0[(*it)->global_index] = (*it)->get_value(time);
            }
        }
    }

    void ExplicitSystem::apply_external_forces(double time) {
        for (ForceList::const_iterator it = external_forces.begin(); it != external_forces.end(); ++it) {
            forces[(*it)->global_index] = (*it)->get_value(time);
        }
    }

} // namespace explicit_fea