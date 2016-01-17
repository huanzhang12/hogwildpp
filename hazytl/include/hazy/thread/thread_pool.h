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

#ifndef HAZY_THREAD_THREAD_POOL_H
#define HAZY_THREAD_THREAD_POOL_H

#include <pthread.h>

#include "hazy/thread/barrier_t.h"

namespace hazy {
namespace thread {

/*! \brief A structure to allow threads to sychornize themselves with the pool.
 */
struct ThreadMeta {
  unsigned thread_id; //!< [0, 1, ..., N]
  bool exit_flag; //!< set to true to have thread exit safely
  barrier_t *ready; //!< wait for a task to work on
  barrier_t *finished; //!< wait after the task is done
  void *tpool; //!< callback into the threadpool that created it
};

/*! \brief A simple homogenous worker thread pool
 */
class ThreadPool {
 public:
  /*! \brief Creates a new thread pool, must call Init() before use.
   * \param n_threads Number of threads
   */
  explicit ThreadPool(unsigned n_threads) : n_threads_(n_threads), 
      threads_(NULL) { }

  virtual ~ThreadPool();

  /*! \brief create the treads for this pool, call exaclty once before use */
  void Init();

  /*! \brief Invoke the function on the task in each thread
   * The thread id and total number of threads are passed to the function.
   * Call Wait() to wait for ALL threads to finish before calling Execute again.
   * Thread ID is zero based index.
   * \param task the task for all threads to share
   * \param hook the function to execute in each thread
   */
  template <class Task>
  void Execute(Task &task, void (*hook)(Task&, unsigned, unsigned));

  /*! \brief block until all threads are done with the task
   * Call exaclty once after each Execute(...)
   */
  void Wait();

  /*! Call to tear down the threads and join them, only call once before d'tor.
   * \note do not call if threads are execting. 
   */
  void Join();

  /*! \brief for private use
   * allows threads to call back in to execute
   */
  inline void ThreadLoop(ThreadMeta &meta);

  /*! \brief for private use
   * Allows the threads to call back into this class
   */
  template <class Task>
  void ThreadCallback(unsigned thread_id);

  /*! \brief Returns the number of threads in this pool.
   * \return number of threads in this pool
   */
  unsigned ThreadCount() { return n_threads_; }

 private:
  unsigned n_threads_; //!< number of thrads in the pool
  barrier_t ready_; //!< wait for a task to work on
  barrier_t finished_; //!< wait when done working on a task
  ThreadMeta *metas_; //!< one for each thread, thread-specific meta data
  pthread_t *threads_; //!< pthead objects, one per thread
  bool ready_flag_; //!< ready for a call to Execute()
  void *arg_; //!< the Task& argument from Execute, cast to void*
  void *hook_; //!< the hook* from Execute, cast to void*
  void (*callback_)(ThreadPool&, unsigned); //!< template nonsense callback

  /*! \brief tell all the threads to exit safely (at some point) */
  void SetExitFlags(bool flag) {
    for (unsigned i = 0; i < n_threads_; i++)
      metas_[i].exit_flag = flag;
  }

  ThreadPool(const ThreadPool&);
  void operator=(const ThreadPool&);
};

} // namespace concur
} // namespace hazy
#endif
