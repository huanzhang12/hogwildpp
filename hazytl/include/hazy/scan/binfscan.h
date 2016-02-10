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

#ifndef HAZY_HOGWILD_FSCAN_BINFSCAN_H
#define HAZY_HOGWILD_FSCAN_BINFSCAN_H

#include <cstdio>
#include <cmath>
#include <algorithm>
#include <inttypes.h>

#include "hazy/types/tuple.h"

namespace hazy {
namespace scan {

//! Scans a file of binary tuples (row, col, value)
/*! The file format is:
 * Bits 0-32: signed integer (endianess unsepcificed) which is the total number
 * of entries which follow.
 * Followed immediately by Entry structs of unspecified padding and endianess.
 *
 * Reads a chunk from the file into temporary buffer; allows iteration one 
 * Entry at a time of the file. Also supports resetting back to the start of 
 * the file. The file is read in the order it appears on disk. 
 *
 * XXX Undefined behavior if file is not formatted properly or if the the file
 * conatins no Entries.
 */
class BinaryFileScanner {
 public:
  //! Opens the file but does not read until first call of Peek() or Next()
  /*! Creates the scanner, ready for calling Next() or Peek(). 
   * \param fname a path to the file on disk to open
   */
  BinaryFileScanner(const char *fname);
  ~BinaryFileScanner();

  //! True if there are more Entries to read by calling Next()
  /*! If !HasNext(), then the behavior of Peek() or Next() is undefined.
   * If HasNext(), then Peek() and Next() will return the next Entry in 
   * the file. Check HasNext() is a constant time operation.
   * \return True if there is more read with Next()
   */
  inline bool HasNext() const;

  //! Returns the next Entry, or dies if !HasNext()
  /*! Returns the next tuple. The file will be iterated in order it written
   * on disk. This is an important and assumed behavior, as the file may be 
   * sorted. Undefined/bad behavior if !HasNext().
   * Do not hold on to the Entry that is returned; immediately make a copy of 
   * it. The value of the returned entry is only valid until the time Peek()
   * or Next() or Reset() is called.
   * \return The next entry in the file, that must be copied to a new location
   * before use.
   */
  inline const types::Entry& Next();

  inline size_t BulkNext(types::Entry *arr, size_t max);
  
  //! Returns the next Entry, or dies if !HasNext()
  /*! Returns the next entry without advancing the iterator. Same rules
   * as Next().
   * \return The next entry in the file, that must be copied to a new location
   * before use.
   */
  inline const types::Entry& Peek();

  //! Returns the maximum value of Entry::col of all entries returned by Next()
  inline int MaxColumn() const { return max_col_; }

  //! Resets the scan
  /*! Resets the iterator to be at the start of the file. Can be called at any
   * time, even if not all the way through the file.
   */
  void Reset();

 private:
  FILE* fh_; //!< Open file handle to the binary file on disk
  const char *fname_; //!< Name of the file bound to fh_
  size_t posn_; //!< the current position in array_ of the Next entry
  size_t size_; //!< the number of entries in the array_
  int max_col_; //!< largest value of Entry::col seen so far
  size_t total_; //!< total remaining Entries to read from the file handle
  size_t buf_size_; //!< total size of array_
  types::Entry* array_; //!< Array of Entries, a chunk of the file. Points into
                      //!< buf_
};

} // namespace scan
} // namespace hazy

#include "binfscan.hxx"
#endif
