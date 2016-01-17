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

#ifndef HAZYHTL_TYPES_TUPLE_H
#define HAZYHTL_TYPES_TUPLE_H

namespace hazy {
namespace types {

/*! \brief An entry in a DIM mode tensor
 */
template <unsigned DIM>
struct Tuple {
  int posn[DIM];
  double rating;

  Tuple() { }

  explicit Tuple(const Tuple &o) {
    for (int i = 0; i < DIM; i++) {
      posn[i] = o.posn[i];
    }
    rating = o.rating;
  }
 private:
  void operator=(const Tuple<DIM> &);
};

/*! \brief A single entry in a matrix
 */
template <>
struct Tuple<2> {
  int posn[]; //!< dubious way to keep consistency with Tupe<DIM>
  int row, col;
  double rating;

  Tuple() { }

  Tuple(const Tuple<2> &o) {
    row = o.row;
    col = o.col;
    rating = o.rating;
  }

  inline void operator=(const Tuple<2> &o) {
    row = o.row;
    col = o.col;
    rating = o.rating;
  }
};

/* \brief An alias for a matrix entry */
typedef Tuple<2> Entry;

} // namespace types
} // namespace hazy
#endif
