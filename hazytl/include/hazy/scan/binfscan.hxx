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


#include <string.h>
#include <assert.h>

hazy::scan::BinaryFileScanner::BinaryFileScanner( const char *fname) 
      : fh_(NULL), fname_(fname), posn_(0), size_(0), max_col_(0), 
      array_(NULL) {
  if ((fh_ = fopen(fname_, "rb")) == NULL) {
    char buf[1024];
    sprintf(buf, "fopen failed for %s", fname_);
    perror(buf);
    exit(-1);
  }
  buf_size_ = 1024 * 1024;
  array_ = new types::Entry[buf_size_];
  Reset();
}

bool hazy::scan::BinaryFileScanner::HasNext() const {
  return total_ > 0;
}

const hazy::types::Entry& 
hazy::scan::BinaryFileScanner::Peek() {
  using hazy::types::Entry;
  assert(HasNext());
  if (posn_ < size_) {
    return array_[posn_];
  }
  // page in from file, this may not fill our buffer completely
  size_t toread = std::min(buf_size_, total_);
  assert(toread > 0);
  size_ = fread(array_, sizeof(Entry), toread, fh_);
  assert(size_ == toread);
  assert(size_ <= total_);
  // clear our buffer
  posn_ = 0;
  return array_[posn_];
}

const hazy::types::Entry& 
hazy::scan::BinaryFileScanner::Next() {
  const hazy::types::Entry &nxt = Peek();
  posn_++;
  total_--;
  max_col_ = std::max(max_col_, nxt.col);
  return nxt;
}

void hazy::scan::BinaryFileScanner::Reset() {
  using hazy::types::Entry;
  fseek(fh_, 0x0, SEEK_SET);
  // FIXME should change fileformat to make it size_t not int
  assert(fread(&total_, sizeof(uint64_t), 1, fh_) == 1);
  size_ = 0;
  posn_ = 0;
}

size_t hazy::scan::BinaryFileScanner::BulkNext(types::Entry *arr, size_t max) {
  Peek();
  size_t load = size_ - posn_;
  if (load > max) {
    load = max;
  }
  memcpy(arr, &array_[posn_], load);
  posn_ += load;
  total_ -= load;
  return load;
}


