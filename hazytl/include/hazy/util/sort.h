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

// Part of the Hazy Project
// Author: Victor Bittorf, bittorf [at] cs.wisc.edu

#ifndef HAZY_UTIL_SORT_H
#define HAZY_UTIL_SORT_H

#include <assert.h>

namespace hazy {
namespace util {

namespace __sort_h {
class SimpleLT {
 public:
  template <class T>
  static inline bool LessThan(T const&a, T const&b) {
    return a < b;
  }
};
}

template <class Comp, typename T>
void QuickSort(T *arr, size_t len) {
  if (len <= 1)
    return;
  //printf("recur: %llu\n", len);

  size_t pivot_index = 0;
  long long gt = len-1;
  T temp;
  T piv = arr[pivot_index];
  // swap pivot to the end
  arr[pivot_index] = arr[gt];
  arr[gt] = piv;
  gt--;
  for (long long i = 0; i <= gt && gt >= 0; ) {
    if (Comp::LessThan(piv, arr[gt])) { gt--; continue; }
    if (Comp::LessThan(piv, arr[i])) {
      // Greater than pivot
      temp = arr[gt];
      arr[gt] = arr[i];
      arr[i] = temp;
      gt--;
    } else {
      i++;
    }
  }
  //assert(Comp::LessThan(piv, arr[gt+1]));
  temp = arr[gt+1];
  arr[gt+1] = piv;
  arr[len-1] = temp;

  // [0,gt] <= piv; gt+1==piv; [gt+1,len-1] >= piv
  assert(gt >= -1);
  long long uniq = gt;
  while (uniq > 0 && arr[uniq] == piv) { uniq--; }
  if (uniq > 0) {
    QuickSort<Comp>(arr, uniq+1);
  }
  if (len - (gt+2) > 1) {
    QuickSort<Comp>(&arr[gt+2], len - (gt+2));
  }
}
template <typename T>
void QuickSort(T *arr, size_t len) {
  return QuickSort<__sort_h::SimpleLT>(arr, len);
}

} // namespace util
} // namespace hazy
#endif
