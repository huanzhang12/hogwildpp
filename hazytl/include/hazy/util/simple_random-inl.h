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

#ifndef HAZY_UTIL_SIMPLE_RANDOM_INL_H
#define HAZY_UTIL_SIMPLE_RANDOM_INL_H

#include <cstdlib>
#include <string.h>

// See for documenation
#include "simple_random.h"
#include "nsec_time.h"

namespace hazy {
namespace util {

void SimpleRandom::SeedByTime() {
  unsigned int seed = (unsigned int) CurrentNSec();
  SetSeed(seed);
}

void SimpleRandom::SetSeed(unsigned int seed) {
  seed_ = seed;
  srand(seed_);
  srand48(seed_);
}

unsigned int SimpleRandom::RandInt(unsigned int nMax) {
  return (unsigned int) (((double) nMax) * (rand_r(&seed_) / (RAND_MAX + 1.0)));
}

double SimpleRandom::RandDouble() { return drand48(); } 

SimpleRandom& SimpleRandom::GetInstance() {
  if (inst_ == NULL) {
    inst_ = new SimpleRandom();
  }
  return *inst_;
}

template <class T>
void SimpleRandom::LazyPODShuffle(T* arr, unsigned int len) {
  char swap[sizeof(T)];
  T *temp = reinterpret_cast<T*>(swap);
  for (unsigned int i = 0; i < len; i++) {
    unsigned int swapi = RandInt(len);
    if (swapi == i) { 
      continue;
    }
    memcpy(temp, &arr[i], sizeof(T));
    memcpy(&arr[i], &arr[swapi], sizeof(T));
    memcpy(&arr[swapi], temp, sizeof(T));
  }
}

} // namespace util
} // namespace hazy
#endif
