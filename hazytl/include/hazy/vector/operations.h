// Copyright 2012 Chris Re, Victor Bittorf
//
 //Licensed under the Apache License, Version 2.0 (the "License");
 //you may not use this file except in compliance with the License.
 //You may obtain a copy of the License at
 //    http://www.apache.org/licenses/LICENSE-2.0
 //Unless required by applicable law or agreed to in writing, software
 //distributed under the License is distributed on an "AS IS" BASIS,
 //WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 //See the License for the specific language governing permissions and
 //limitations under the License.

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_VECTOR_VECTOROPS_H
#define HAZY_VECTOR_VECTOROPS_H

#include <cstdlib>
#include <cmath>
#include <algorithm>

#include "fvector.h"
#include "svector.h"

namespace hazy {
namespace vector {

/*! \brief Projects the given SVector using the index mask into the out array
 * \param s the vector to project
 * \param indexes the indicies to write
 * \param out the array to write to
 */
template <typename T, typename int_T>
void inline Project(SVector<T> const&s, FVector<int_T> const &indexes, T *out);

/*! \brief Preforms a simplex projection on the given vector
 * \param v vector to project (and modify)
 */
template <typename float_t>
void inline SimplexProject(FVector<float_t> &v);

/*! \brief sets all values of the vector to zero
 * \param v vector to zero
 */
template <typename float_t>
void inline Zero(FVector<float_t> &v);

/*! \brief Sets any element that is less than zero to be zero
 * \param v vector to threshold
 */
template <typename float_t>
void inline ThresholdZero(SVector<float_t> &v);

/*! \brief Copy one vector into another vector
 * Call must ensure that the backing memory of out is big enough
 * \param u vector to copy from
 * \param out vector to copy into
 */
template <typename float_u>
void inline CopyInto(FVector<float_u> const &u, FVector<float_u> &out);

} // namespace vector
} // namespace hazy
#endif
