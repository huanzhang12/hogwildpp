
#ifndef HOGWILD_TEST_FAKE_SCAN_INL_H
#define HOGWILD_TEST_FAKE_SCAN_INL_H

#include "hazy/vector/fvector.h"
#include "hazy/hogwild/hogwild_task.h"

namespace hazy {
namespace hogwild {

template <class T>
class FakeScan {
 public:
  FakeScan(int passes, int block_size) : done_(0), passes_(passes),
      block_size_(block_size) {
    fake.ex.size = block_size;
    fake.perm.size = block_size;
  }

  bool HasNext() { return done_ < passes_; }

  ExampleBlock<T>& Next() {
    done_++;
    return fake;
  }

  void Reset() { done_ = 0; }

  int done_;
  int passes_;
  int block_size_;
  ExampleBlock<T> fake;
};

} // namespace hogwild
} // namespace hazy
#endif
