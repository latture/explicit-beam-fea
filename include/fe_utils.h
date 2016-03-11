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

#ifndef EXPLICIT_FEA_FE_UTILS_H
#define EXPLICIT_FEA_FE_UTILS_H

#include "beam_element.h"
#include <memory>
#include <type_traits>

namespace explicit_fea {
    /**
     * @brief Updates rotation matrices to the current element's geometry based on the
     * associated vector of `Node`s.
     * 
     * @param nodes[in]               List of nodal coordinates.
     * @param elem[in]                Current element.
     * @param rotation[in/out]            Square matrix that transforms elemental matrices from local to global coordinates.
     * @param rotation_transposed[in/out] Transpose of the aforementioned rotation matrix.
     */
    void update_rotation(const std::vector<Node>& nodes,
                         const BeamElement& elem,
                         LocalMatrix& rotation,
                         LocalMatrix& rotation_transposed);

    /**
     * Estimates the maximum stable time step.
     * @param  nodes[in] List of nodal coordinates.
     * @param  elems[in] List of finite elements associated with the nodal coordinates.
     * @return           Stable time step.
     */
    double estimate_stable_timestep(const std::vector<Node>& nodes, const std::vector<std::unique_ptr<BeamElement>>& elems);

    /**
     * Object that can be used to either integral or floating point values.
     */
    template <typename T>
    class ValueCompare {
    public:
        /**
         * @brief Constructor
         * @param epsilon Allowed absolute tolerance when comparing 2 values.
         */
        ValueCompare(T epsilon=1.e-14) : epsilon(epsilon) {}

        /**
         * @brief Returns whether input `a` is less than input `b`
         * @details If `a` and `b` are not integral types then the
         * values are compared with an absolute tolerance given by
         * the member variable `epsilon`
         *
         * @param[in] a First value
         * @param[in] b Second value
         * @return `a < b`
         */
        bool lessThan(const T &a, const T &b) const {
            if (std::is_integral<T>::value) {
                return a < b;
            }
            else
                return (b - a) > ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * epsilon);
        }

        /**
         * @brief Returns whether input `a` is greater than input `b`
         * @details If `a` and `b` are not integral types then the
         * values are compared with an absolute tolerance given by
         * the member variable `epsilon`
         *
         * @param[in] a First value
         * @param[in] b Second value
         * @return `a > b`
         */
        bool greaterThan(const T &a, const T &b) const {
            if (std::is_integral<T>::value) {
                return a > b;
            }
            else
                return (a - b) > ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * epsilon);
        }

        /**
         * @brief Returns whether input `a` is equal to input `b`
         * @details If `a` and `b` are not integral types then the
         * values are compared with an absolute tolerance given by
         * the member variable `epsilon`
         *
         * @param[in] a First value
         * @param[in] b Second value
         * @return `a == b`
         */
        bool equal(const T &a, const T &b) const {
            if (std::is_integral<T>::value) {
                return a == b;
            }
            else
                return std::abs(a - b) <= ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * epsilon);
        }
    private:
        const T epsilon;
    };
} // namespace explicit_fea

#endif //EXPLICIT_FEA_FE_UTILS_H
