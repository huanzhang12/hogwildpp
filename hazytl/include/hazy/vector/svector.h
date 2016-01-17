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

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_TYPES_SVECTOR_H
#define HAZY_TYPES_SVECTOR_H
#ifndef nullptr
#define nullptr NULL
#endif

#include <cstdlib>
#include <inttypes.h>

namespace hazy {
namespace vector {

/*! \brief Sparse vector which contains implied zero entries
 * A representation of a vector where elements are encoded as a pair of
 * (index, value). An SVector stores exactly size such pairs, where any pair
 * pair not explicitly stored is assumed to be zero. Pair i is encoded as
 * (SVecetor::index[i], SVector::values[i]) for all 0<= i < size;
 * All memory is managed by callers.
 */
template <typename T>
struct SVector {
 public:
  uint64_t size; //!< number of pairs, length of index & values
  int *index; //!< accending array of indicies
  T *values; //!< array of values for the corresponding indicies

  //! Create an SVector of the given size, backed by the provided arrays
  /*! Creates a new SVector backed by the given arrays, of the given size.
   * Caller must manage the memory.
   * \param vals the values of the vector, array of length s
   * \param idx the indicies of the given values, array of length s
   * \param s the size of the vector
   */
  explicit SVector(T * vals, int * idx, uint64_t s) :
      size(s), index(idx), values(vals) { }

  //! Create an SVector of the given size, backed by the provided arrays
  /*! Creates a new SVector backed by the given arrays, of the given size.
   * Caller must manage the memory.
   * \param vals the values of the vector, array of length s
   * \param idx the indicies of the given values, array of length s
   * \param s the size of the vector
   */
  explicit SVector(uint64_t s, T * vals, int * idx) :
      size(s), index(idx), values(vals) { }
  
  /*! \brief Create a new SVector copy of an existing one
   * \param other the vector to copy
   */
  SVector(SVector<T> const & other) {
    size = other.size;
    index = other.index;
    values = other.values;
  }

  /*! \brief Create an empty vector. */
  SVector() : size(0), index(nullptr), values(nullptr) { }

  /*! \brief Sets this vector to a copy of the given vector
   * \param o vector to copy
   * \return this vector
   */
  SVector<T>& operator=(SVector<T> const &o) {
    size = o.size;
    index = o.index;
    values = o.values;
    return *this;
  }
};

} // namespace vector
} // namespace hazy
#endif
