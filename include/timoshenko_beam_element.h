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

#ifndef EXPLICIT_BEAM_FEA_TIMOSHENKO_BEAM_ELEMENT_H
#define EXPLICIT_BEAM_FEA_TIMOSHENKO_BEAM_ELEMENT_H

#include "beam_element.h"

namespace explicit_fea {

    /**
     * @brief Implementation of the Timoshenko beam element
     * @details The Timoshenko beam element is intended for shear-deformable beams and 
     * operates under that assumption the plan cross-sections remain plane but can rotate relative
     * to the neutral axis during deformation.
     */
    class TimoshenkoBeamElement : public BeamElement {
    public:
        /**
         * @brief Constructor
         * @details constructs a beam element of specified properties and refencing the 
         * specifed node numbers of an associated node list.
         *
         * @param[in] nn1 Index of the first node of the element in an associated node list (see `NodeList`).
         * @param[in] nn2 Index of the second node of the element in an associated node list (see `NodeList`).
         * @param[in] props Elemental properties.
         */
        TimoshenkoBeamElement(int nn1, int nn2, const Props& props);

        /**
         * @brief Returns the elemental stiffness matrix.
         * @details The elemental stiffness matrix is a square matrix of size `2N`x`2N`.
         * Similar formulation to `EulerBernoulliBeamElement::calculate_stiffness_matrix`
         * with an additional "corrective" term to allow for shear deformation.
         *
         * @param[in] nodes List of nodal coordinates, where each node cooresponds to a given \f$ (x, y, z) \f$ point.
         * @return Elemental stiffness matrix.
         */
        LocalMatrix calculate_stiffness_matrix(const std::vector<Node>& nodes);
        
        /**
         * @brief Returns the inverse mass matrix.
         * @details The consistent mass matrix is calculated, derived from the same shape functions used with calculating
         * the elemental stiffness matrix.
         * 
         * @param[in] nodes List of nodal coordinates, where each node cooresponds to a given \f$ (x, y, z) \f$ point.
         * @return Elemental inverse mass matrix.
         */   
        LocalMatrix calculate_inv_mass_matrix(const std::vector<Node>& nodes);
    };

} // namespace explicit_fea

#endif //EXPLICITBEAMFEA_TIMOSHENKO_BEAM_ELEMENT_H
