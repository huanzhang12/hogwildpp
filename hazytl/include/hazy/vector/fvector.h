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

#ifndef HAZY_TYPES_FVECTOR_H
#define HAZY_TYPES_FVECTOR_H
#ifndef nullptr
#define nullptr NULL
#endif

#include <cstdlib>
#include <inttypes.h>

namespace hazy {
namespace vector {

/*! \brief This is a container for a size and pointer of a block of memory.
 * Memory must be managed by caller. 
 */
template <typename T>
class FVector {
 public:
  uint64_t size; //!< the length of the vector
  T *values; //!< pointer to the backing array

  //! Create a new, empty vector
  FVector() : size(0), values(nullptr) { }

  /*! \brief Create a vector of length s and of elements at p.
   * Create a new vector. Caller must manage the memory.
   * \param p the values of the vector, must be of (atleast) length s
   * \param s the length of the vector
   */
  explicit FVector(T *p, uint64_t s) : size(s), values(p) { }

  /*! \brief Create a vector of length s and of elements at p.
   * Create a new vector. Caller must manage the memory.
   * \param s the length of the vector
   * \param p the values of the vector, must be of (atleast) length s
   */
  explicit FVector(uint64_t len, T *p) : size(len), values(p) { }

  /*! \brief Use an existing FVector ot initialize this one
   * \param o existing vector to copy
   */
  FVector(FVector<T> const &o) {
    size = o.size;
    values = o.values;
  }

  /*! \brief Sets this vector to be a copy of the other vector
   * \param o the vector to copy
   */
  FVector<T>& operator=(FVector<T> const& o) {
    size = o.size;
    values = o.values;
    return *this;
  }

  inline T& operator[](size_t idx) {
    return values[idx];
  }
};

} // namespace vector
} // namespace hazy
#endif
