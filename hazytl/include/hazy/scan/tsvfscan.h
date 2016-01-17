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

// Hogwild!, part of the Hazy Project
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)
// Original Hogwild! Author: Chris Re (chrisre [at] cs.wisc.edu)             
#ifndef HAZY_HOGWILD_FSCAN_TSVFSCAN_H
#define HAZY_HOGWILD_FSCAN_TSVFSCAN_H
#include <cmath>

#include <cstdio>
#include <iostream>
#include <fstream>

#include "hazy/types/tuple.h"

namespace hazy {
namespace scan {

//! Enumerates tuples from a TSV file
template <int OFFSET, int DIM>
class OffsetTSVFileScanner {
 public:
  //! Creates a new scanner for the file
  explicit OffsetTSVFileScanner(const char *fname);

  //! True if Next() will return a valid value
  inline bool HasNext();
  //! Returns the next Entry, or dies if !HasNext()
  inline const types::Tuple<DIM>& Next();
  //! Returns the next Entry, or dies if !HasNext()
  inline const types::Tuple<DIM>& Peek();
  //! Returns the max column seen so far by Next()
  inline int MaxColumn() const { return max_col_; }
  //! Resets the scan
  void Reset();

 private:
  FILE* fh;
  const char *fname_; //!< memo of the file name
  bool peekable_; //!< flag set to true if entry_ has valid data
  types::Tuple<DIM> entry_; //!< the item to return next
  int max_col_; //!< biggest column seen so far
};

//! Zero-offset TSV File Scanner
typedef OffsetTSVFileScanner<0, 2> TSVFileScanner;
//! TSV Scanner offsets rows & columns by -1, to match matlab
typedef OffsetTSVFileScanner<-1, 2> MatlabTSVFileScanner;

} // namespace scan
} // namespace hazy

#include "tsvfscan.hxx"
#endif
