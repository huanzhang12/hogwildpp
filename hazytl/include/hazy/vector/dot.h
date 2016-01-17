// Copyright 2012 Chris Re, Victor Bittorf
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_VECTOR_DOT_H
#define HAZY_VECTOR_DOT_H

namespace hazy {
namespace vector {

/*! \brief Compute the dot product, missing entries in the SVector are assumed 0.
 * \return the dot product.
 */
template <typename float_u, typename float_v>
float_u inline Dot(FVector<float_u> const&, SVector<float_v> const&);


/*! \brief Compute the dot product, assumes the vectors are the same length.
 * \return the dot product.
 */
template <typename float_u, typename float_v>
float_u inline Dot(FVector<float_u> const&, FVector<float_v> const&);

/*! \brief Compute the dot product, assumes missing entries are zero.
 * \return the dot product.
 */
template <typename float_u, typename float_v>
float_u inline Dot(SVector<float_u> const&, SVector<float_v> const&);

} // namespace vector
} // namespace hazy
#endif
