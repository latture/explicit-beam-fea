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

#ifndef EXPLICIT_FEA_PRESCRIBED_VALUE_H
#define EXPLICIT_FEA_PRESCRIBED_VALUE_H

#include "containers.h"

namespace explicit_fea {

    /**
     * @brief Nodal value set by the user.
     * @details Values can be prescribed via external forces (`explicit_fea::Force`)
     * or boundary conditions (`explicit_fea::BC`). All prescribed values act on a
     * nodal degree of freedom. This class provides a common structure
     * for specifying which global degree of freedom is prescribed.
     */
    struct PrescribedValue {
        /**
         * @brief Returns the prescribed nodal value at the specified time.
         * @details Override this method when creating a valid prescribed value.
         * @sa `explicit_fea::BC`, `explicit_fea::Force`
         */
        virtual double get_value(double time) = 0;

        /**
         * The index of the node to constrain
         */
        const int node;
        /**
         * The index of the dof to constrain. The `explicit_fea::DOF` enum can be used for specification or the integer
         * values can be used directly `0==d_x`, `1==d_y`, ...
         *
         * @sa `explicit_fea::DOF`
         */
        const int dof;

        /**
         * @brief Index of the prescribed value in global context.
         * @details The global index associated with a prescribed nodal value can be found from
         * `exlicit_fea::DOF::NUM_DOFS * node + dof`, i.e. the number of degrees of freedom per node
         * mulitplied by the node number plus the local DOF being constrained.
         */
        const int global_index;

    protected:
        /**
         * @brief Constructor
         * @param[in] node `int`. The index of the node.
         * @param[in] dof `int`. Degree of freedom to constrain (See `explicit_fea::DOF`).
         */
        PrescribedValue(int _node, int _dof) : node(_node), dof(_dof), global_index(DOF::NUM_DOFS * _node + _dof) {}
    };
}

#endif //EXPLICIT_FEA_PRESCRIBED_VALUE_H
