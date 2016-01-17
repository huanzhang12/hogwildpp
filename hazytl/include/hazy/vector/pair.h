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

#ifndef HAZY_VECTOR_PAIR_H
#define HAZY_VECTOR_PAIR_H

#include <cstdlib>
#include <inttypes.h>

#include "fvector.h"

namespace hazy {
namespace vector {

template <typename T>
struct Pair {
  T value;
  int index;

  Pair() { }

  Pair(T v, int i) : value(v), index(i) { }

  Pair(Pair const& p) {
    value = p.value;
    index = p.index;
  }

  inline Pair& operator=(Pair const& p) {
    value = p.value;
    index = p.index;
  }
};

typedef FVector<Pair<double> > PVectorDouble;

} // namespace vector
} // namespace hazy

#endif
