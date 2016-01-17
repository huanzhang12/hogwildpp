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

#ifndef HAZY_VECTOR_DOT_INL_H
#define HAZY_VECTOR_DOT_INL_H

#include <cmath>

#include "hazy/vector/dot.h"

// See hazy/vector/dot.h for documentation

namespace hazy {
namespace vector {

template <typename T, typename Vec>
void inline MatrixVectorMultiply(FVector<SVector<T> > const &mat, 
                                 Vec const &vec,
                                 FVector<T> &out) {
  // Note: Assume mat is padded with zeros
  out.size = mat.size;
  for (size_t i = 0; i < mat.size; i++) {
    out.values[i] = Dot(vec, mat.values[i]);
  }
}

template <typename float_u, typename float_v>
float_u inline Dot(FVector<float_u> const& u, SVector<float_v> const& v) {
  float_u p = 0.0;
  float_u const * const /*__restrict__*/ uvals = u.values;
  float_v const * const /*__restrict__*/ vvals = v.values;
  int const * const /*__restrict__*/ vidx = v.index;
  for (size_t i = v.size; i-- > 0; ) {
    p += uvals[vidx[i]] * vvals[i]; 
  }
  return p;
}

template <typename float_u, typename float_v>
float_u inline Dot(FVector<float_u> const &u, FVector<float_v> const &v) {
  float_u p = 0.0;
  float_u const * const /*__restrict__*/ uvals = u.values;
  float_v const * const /*__restrict__*/ vvals = v.values;
  for (size_t i = v.size; i-- > 0; ) {
    p += uvals[i] * vvals[i]; 
  }
  return p;
}

template <typename float_u, typename float_v>
float_u inline Dot(SVector<float_u> const& x, SVector<float_v> const& y) {
  long xi = (long)x.size - 1;
  long yi = (long)y.size - 1;
  float_u ret = 0.0;
  while(xi >= 0 && yi >= 0) {
    int diff = x.index[xi] - static_cast<int>(y.index[yi]);
    if(diff == 0) {
      ret += x.values[xi]*y.values[yi];
      xi--; yi--;
    } else {  
      if(diff > 0) {xi--;}  else { yi--; }
    }
  }
  return ret;
}

} // namespace vector
} // namespace hazy
#endif
