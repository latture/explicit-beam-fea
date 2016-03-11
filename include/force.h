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

#ifndef EXPLICITBEAMFEA_FORCE_H
#define EXPLICITBEAMFEA_FORCE_H

#include <memory>
#include "prescribed_value.h"

namespace explicit_fea {
    /**
     * @brief Parent class of external forces applied to nodal degrees of freedom.
     * @details Inherit from `Force` and implement the `double get_value(double time)` method to create a Force.
     * The `get_value` method allows the user to have external forces that are a function of time.
     *
     * @sa `explicit_fea::ConstantForce`
     */
    struct Force : public PrescribedValue {

        /**
         * @brief Constructor
         * @param[in] node `int`. The index of the node.
         * @param[in] dof `int`. Degree of freedom to constrain (See `explicit_fea::DOF`).
         */
        Force(int _node, int _dof) : PrescribedValue(_node, _dof) { };
    };

    typedef std::vector<std::unique_ptr<Force>> ForceList;

    struct ConstantForce : public Force {
        /**
         * @brief Constructor
         * @param[in] node `int`. The index of the node.
         * @param[in] dof `int`. Degree of freedom to constrain (See fea::DOF).
         * @param[in] value `double`. The prescribed value for the force.
         */
        ConstantForce(int _node, int _dof, double _value) : Force(_node, _dof), value(_value) {};

        double get_value(double time) {
            return value;
        }

    private:
        const double value;
    };

} // namespace explicit_fea

#endif //EXPLICITBEAMFEA_FORCE_H
