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

#ifndef EXPLICIT_FEA_BCS_H
#define EXPLICIT_FEA_BCS_H

#include <memory>
#include "prescribed_value.h"

namespace explicit_fea {
    /**
     * @brief Parent boundary condition class.
     * @details Inherit from `BC` and implement the `double get_value(double time)` method to create a BC.
     * The `get_value` method allows the user to have boundary conditions that are a function of time.
     *
     * @sa `ConstantBC`
     */
    struct BC : public PrescribedValue {
        /**
         * Enum used to specify whether the boundary condition applies to
         * the velocity or displacement of the nodal degree of freedom.
         */
        enum Type {
            DISPLACEMENT,/**<Apply boundary condition to displacement of nodal degree of freedom.*/
            VELOCITY /**<Apply boundary condition to velocity of nodal degree of freedom.*/
        };

        /**
         * @brief Constructor
         * @param[in] node The index of the node.
         * @param[in] dof Degree of freedom to constrain (See `explicit_fea::DOF`).
         * @param[in] type The type of boundary condition, can be either `BC::VELOCITY` or `BC::DISPLACEMENT`.
         */
        BC(int _node, int _dof, Type _type)
                : PrescribedValue(_node, _dof), type(_type) { };

        virtual ~BC() {};

        const Type type;/**<Type of boundary condition. Can be either velocity of displacement.*/
    };

    typedef std::vector<std::unique_ptr<BC>> BCList;

    /**
     * @brief Boundary condition that is constant in time.
     * @details Calling `get_value` will always return the initial `value` specified.
     */
    struct ConstantBC : public BC {
        /**
         * @brief Constructor
         * @param[in] node The index of the node.
         * @param[in] dof Degree of freedom to constrain (See `explicit_fea::DOF`).
         * @param[in] value The prescribed value for the boundary condition.
         * @param[in] type The type of boundary condition, can be either `BC::VELOCITY` or `BC::DISPLACEMENT`.
         */
        ConstantBC(int _node, int _dof, double _value, Type _type) : BC(_node, _dof, _type), value(_value) {}

        double get_value(double time) {
            return value;
        }

    private:
        const double value;/**<The value to hold the dof at.*/
    };
} // namespace explicit_fea

#endif //EXPLICIT_FEA_BCS_H
