// Copyright 2012 Victor Bittorf, Chris Re
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

// Hazy Template Library 
// Created: October 2012
// Primary Author: Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_UTIL_SIMPLE_RANDOM_H
#define HAZY_UTIL_SIMPLE_RANDOM_H

#ifndef HAZY_UTIL_SIMPLE_RANDOM_INL_H
#error Include -inl.h
#endif

#include <time.h>

namespace hazy {
namespace util {

/*! \brief Wraper for random number generation
 */
class SimpleRandom {
 public:
  /*! \brief Sets the seed all newly constructed RNGs will use by default 
   * \param seed the seed for SimpleRandom objects to use as a default
   */
  static void SetSeed(unsigned int seed);

  /*! \brief Uses the current time as the argument to SetSeed(...) */
  static void SeedByTime();

  /*! \brief Gets the singleton RNG */
  inline static SimpleRandom& GetInstance();
  
  /*! \brief Get a random integer inclusive in [0, nMax)
   */
  inline unsigned RandInt(unsigned int nMax); 

  /*! \brief Get a random double inclusive in [0.0, 1.0]
   */
  inline double RandDouble();

  /*! \brief Puesdo-Shuffles the given array inplace; only handles integer len
   * This is not a rigorious shuffle, it is not rigorious shuffle. This shuffle
   * is designed to be fast and in-place but may not yield all possible
   * orderings or may not uniformly shuffle.
   * \param arr the array of POD types to shuffle
   * \param len the number of elements in arr
   */
  template <class T>
  inline void LazyPODShuffle(T* arr, unsigned int len);

 private:
  /*! \brief last random seed used */
  static unsigned int seed_;
  /*! \brief Singleton instance */
  static SimpleRandom *inst_;

  /*! \brief Generate a new RNG seeded by the default seed set by SetSet(...)
   */
  SimpleRandom() { };
};

SimpleRandom *SimpleRandom::inst_ = NULL;
unsigned int SimpleRandom::seed_ = 0;

} // namespace util
} // namespace hazy
#endif
