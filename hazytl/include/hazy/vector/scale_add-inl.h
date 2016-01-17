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

#ifndef HAZY_VECTOR_SCALE_ADD_INL_H
#define HAZY_VECTOR_SCALE_ADD_INL_H

// See for documentation
#include "hazy/vector/scale_add.h"

namespace hazy {
namespace vector {

template <typename float_u, typename float_v, typename float_>
void inline DropScaleAndAdd(SVector<float_u> &x, SVector<float_v> const& y, 
                   float_ scal) {
  size_t xi = x.size - 1;
  size_t yi = y.size - 1;
  while(xi > 0 && yi > 0) {
    int diff = x.index[xi] - static_cast<int>(y.index[yi]);
    if(diff == 0) {
      x.values[xi] += y.values[yi] * scal;
      xi--; yi--;
    } else {
      if(diff > 0) {xi--;}  else { yi--; }
    }
  }
}

template <typename float_u, typename float_v, typename float_>
void inline ScaleAndAdd(FVector<float_u> &u, SVector<float_v> const &v,
                                             float_ const& s) {
  float_u * const __restrict__ uvals = u.values;
  float_v const * const __restrict__ vvals = v.values;
  int const * const __restrict__ vidx = v.index;
  for (size_t i = v.size; i-- > 0; ) {
    uvals[vidx[i]] = uvals[vidx[i]] + vvals[i] * s;
  }
}

template <typename float_u, typename float_v, typename float_>
void inline ScaleAndAdd(FVector<float_u> &u, FVector<float_v> const &v,
                                             float_ const& s) {
  float_u * const __restrict__ uvals = u.values;
  float_v const * const __restrict__ vvals = v.values;
  for (size_t i = v.size; i-- > 0; ) {
    uvals[i] += vvals[i] * s;
  }
}

template <typename float_u, typename float_v, typename float_,
          typename t_float>
void inline ScaleAndAddInto(FVector<float_u> const &u, 
                        FVector<float_v> const &v, float_ const& s,
                        FVector<t_float> &copyout) {
  t_float * const __restrict__ out = copyout.values;
  float_v const * const __restrict__ vvals = v.values;
  float_u const * const __restrict__ uvals = u.values;
  for (size_t i = v.size; i-- > 0; ) {
    out[i] = uvals + vvals * s;
  }
}

template <typename float_u, typename float_v>
void inline ScaleInto(FVector<float_u> const &u, float_v const &s, 
                      FVector<float_u> &copyto) {
  float_u * const __restrict__ arr = copyto.values;
  float_u const * const __restrict__ uarr = u.values;
  for (size_t i = u.size; i -- > 0; ) {
    arr[i] = uarr[i] * s;
  }
  copyto.size = u.size;
}

template <typename float_u, typename float_v>
void inline Scale(FVector<float_u> &u, float_v const &s) {
  float_u * __restrict__ const uarr = u.values;
  for (size_t i = u.size; i-- > 0; ) {
    uarr[i] *= s;
  }
}

} // namespace vector
} // namespace hazy
#endif
