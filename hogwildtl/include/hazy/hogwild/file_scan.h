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

#ifndef HAZY_HOGWILD_HOGWILD_SCAN_H
#define HAZY_HOGWILD_HOGWILD_SCAN_H

#include "hazy/util/simple_random-inl.h"

#include "hazy/vector/fvector.h"
#include "hazy/thread/thread_pool-inl.h"
#include "hazy/hogwild/hogwild_task.h"

namespace hazy {
namespace hogwild {

/*! \brief Scan a file using shadow buffering to improve performance
 * \tparam A scanner which produces Example
 * \tparam Example the type of example produced
 */
template <class Scan, class Example>
class FileScan {
 public:
  /*! \brief Create a shadow buffered scan using at most max_mem bytes in memory
   * \param scan A scanner over the file (e.g. BinaryFileScanner)
   * \param max_mem The maximum amount of memory to use in bytes, this is
   *    split between the main buffer and the shadow buffer.
   */
  FileScan(Scan &scan, size_t max_mem) : scan_(scan), has_next_(false), 
    tpool_(1), pageno_(0), max_mem_(max_mem) { }

  ~FileScan() {
    delete [] blk_.ex.values;
    delete [] blk_.perm.values;
    delete [] shadow_blk_.ex.values;
    delete [] shadow_blk_.perm.values;
  }

  /*! \brief Initialize this file scanner. Call exactly once before use.
   */
  void Init() {
    tpool_.Init();
    size_t bufsize = (max_mem_ / (2 * (sizeof(Example) + sizeof(size_t))) + 1);
    blk_.ex.values = new Example[bufsize];
    blk_.perm.values = new size_t[bufsize];
    shadow_blk_.ex.values = new Example[bufsize];
    shadow_blk_.perm.values = new size_t[bufsize];
    task_.max_size = bufsize;
    task_.scan = &scan_;

    Reset();
  }

  /*! \brief Tears down the file scanner, call exactly once after use 
   * but before deconstruction. 
   */
  void Destroy() {
    if (has_next_) {
      tpool_.Wait();
    }
    tpool_.Join();
  }
  
  /*! \brief returns True if there are additoinal pages in this scan.
   * \return true if it is valid to call Next() again
   */
  bool HasNext() { return has_next_; }

  /*! \brief Returns the next block of examples with a permuation.
   * Only call if this HasNext(), undefined otherwise.
   */
  ExampleBlock<Example>& Next() {
    // wait for shadow thread to join
    tpool_.Wait();

    pageno_++;

    has_next_ = false;
    if (scan_.HasNext()) {
      has_next_ = true;
      StartShadowTask();
    }

    return *GetCurrent();
  }

  /*! \brief Reset this scanner to start from the begining of the file.
   * 
   */
  void Reset() { 
    pageno_ = 0;
    if (has_next_) {
      tpool_.Wait();
    }
    scan_.Reset();
    has_next_ = scan_.HasNext();
    if (has_next_) {
      StartShadowTask();
    }
  }

  void StartShadowTask() {
    task_.blk = GetShadow();
    tpool_.Execute(task_, FileScan<Scan, Example>::Delegate);
  }

  ExampleBlock<Example>* GetCurrent() {
    if (pageno_ % 2 == 0) {
      return &blk_;
    }
    return &shadow_blk_;
  }


  ExampleBlock<Example>* GetShadow() {
    if (pageno_ % 2 == 1) {
      return &blk_;
    }
    return &shadow_blk_;
  }
 private:

  struct ShadowTask {
    ExampleBlock<Example> *blk;
    Scan *scan;
    size_t max_size;
  };

  Scan &scan_;
  ExampleBlock<Example> blk_;
  ExampleBlock<Example> shadow_blk_;
  bool has_next_;
  ShadowTask task_;
  hazy::thread::ThreadPool tpool_;
  size_t pageno_;
  size_t max_mem_;

  static void Delegate(ShadowTask &task, unsigned tid, unsigned tot) {
    size_t loaded = 0;
    Scan &scan = *task.scan;
    Example *ex = task.blk->ex.values;
    while ( (scan.HasNext()) && (loaded < task.max_size) ) {
      ex[loaded++] = scan.Next();
    }
    task.blk->ex.size = loaded;
    for (size_t i = 0; i < loaded; i++) {
      task.blk->perm.values[i] = i;
    }
    task.blk->perm.size = loaded;
    util::SimpleRandom &rand = util::SimpleRandom::GetInstance();
    rand.LazyPODShuffle(task.blk->perm.values, loaded);
  }
};

} // namespace hogwild
} // namespace hazy

#endif

