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

#include <assert.h>

namespace hazy {
namespace scan {

template <int OFFSET, int DIM>
OffsetTSVFileScanner<OFFSET, DIM>::OffsetTSVFileScanner(
    const char *fname) : fname_(fname), peekable_(false), max_col_(0) {
  fh = fopen(fname, "r");
  if (fh == NULL) {
    perror("Could not open file");
    assert(false);
  }
  Reset();
}

template <int OFFSET, int DIM>
bool OffsetTSVFileScanner<OFFSET, DIM>::HasNext() {
  if (peekable_) {
    return true;
  }
  return !feof(fh);
}

template <int OFFSET, int DIM>
const hazy::types::Tuple<DIM>& 
OffsetTSVFileScanner<OFFSET, DIM>::Peek() {
  if (peekable_) {
    return entry_;
  }

  if (DIM == 2) {
    assert(fscanf(fh, "%d\t%d\t%lf", &entry_.posn[0], &entry_.posn[1], &entry_.rating) == 3);
  } else {
    for (int i = 0; i < DIM; i++) {
      assert(fscanf(fh, "%d", &entry_.posn[i]) == 1);
    }
    assert(fscanf(fh, "%lf", &entry_.rating) == 1);
  }

  char temp[8];
  while (fscanf(fh, "%1[ \n]s", temp)) { if (feof(fh)) break; }

  peekable_ = true;
  for (int i = 0; i < DIM; i++) {
    entry_.posn[i] += OFFSET;
  }
  return entry_;
}

template <int OFFSET, int DIM>
const hazy::types::Tuple<DIM>& 
OffsetTSVFileScanner<OFFSET, DIM>::Next() {
  Peek(); // make sure entry_ is set
  peekable_ = false;
  return entry_;
}

template <int OFFSET, int DIM>
void OffsetTSVFileScanner<OFFSET, DIM>::Reset() {
  fclose(fh);
  fh = fopen(fname_, "r");
  peekable_ = false;
}

} // namespace fscan
} // namespace hazy
