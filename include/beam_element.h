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

#ifndef EXPLICIT_FEA_ELEMENT_H
#define EXPLICIT_FEA_ELEMENT_H

#include "containers.h"

namespace explicit_fea {

    /**
     * @brief Parent class of Euler-Bernoulli and Timoshenko beam elements.
     * @details Each beam element is comprised of 2 node--these refer to
     * indices in a vector of `Node`'s--and a set of properties. Specific
     * to each element is the method of obtaining the stiffness and
     * consistent mass matrices. These methods must be implemented in
     * the respective derived class.
     *
     * @sa `EulerBernoulliBeamElement`, TimoshenkoBeamElement`
     */
    class BeamElement {
    public:
        virtual ~BeamElement() {};

        /**
         * @brief Returns the elemental stiffness matrix.
         * @details The elemental stiffness matrix should a square matrix of size `2N`x`2N`
         * derived from the shape functions and elemental properties.
         *
         * @param nodes[in] List of nodal coordinates, where each node cooresponds to a given \f$ (x, y, z) \f$ point.
         * @return Elemental stiffness matrix.
         *
         * @sa `EulerBernoulliBeamElement::calculate_stiffness_matrix`, `TimoshenkoBeamElement::calculate_stiffness_matrix`
         */
        virtual LocalMatrix calculate_stiffness_matrix(const std::vector<Node>& nodes) = 0;

        /**
         * @brief Returns the inverse mass matrix.
         * 
         * @param[in] nodes List of nodal coordinates, where each node cooresponds to a given \f$ (x, y, z) \f$ point.
         * @return Elemental inverse mass matrix.
         *
         * @sa `EulerBernoulliBeamElement::calculate_inv_mass_matrix`, `TimoshenkoBeamElement::calculate_inv_mass_matrix`
         */
        virtual LocalMatrix calculate_inv_mass_matrix(const std::vector<Node>& nodes) = 0;

        /**
         * @brief Returns the indices of the nodal coordinates referenced by the beam element.
         * @return Vector of 2 integers corresponding to the nodal indices of the element.
         */
        const Eigen::Vector2i& get_node_numbers() const { return node_numbers; };
        /**
         * @brief Returns the elemental properties
         * @return Elemental properties.
         * @sa Props
         */
        const Props& get_props() const { return props; };

    protected:
        /**
         * Constructor
         *
         * @param[in] nn1 Index of the first node of the element in an associated node list (see `NodeList`).
         * @param[in] nn2 Index of the second node of the element in an associated node list (see `NodeList`).
         * @param[in] props Elemental properties.
         */
        BeamElement(int nn1, int nn2, const Props& props);

        const Props props;/**<Elemental properties.*/
        Eigen::Vector2i node_numbers;/**<Nodal indices refenced by the beam element. These indices are applied to an associated node list.*/
    };

} // namespace explicit_fea

#endif //EXPLICIT_FEA_ELEMENT_H
