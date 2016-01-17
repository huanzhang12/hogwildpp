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


#ifndef HAZY_HOGWILD_HOGWILD2_H
#define HAZY_HOGWILD_HOGWILD2_H

#include "hazy/vector/fvector.h"
#include "hazy/thread/thread_pool-inl.h"

namespace hazy {
namespace hogwild {

/*! \brief Hogwild! parallel executor
 * This is a new version of Hogwild! which allows for more sophiscated
 * aggregations and control of execution. 
 * For example usage, see hazy/hogwild/hogwildexper.h
 */
template <class Model, class Params>
class Hogwild {
 public:
  /*! \brief Sets up for the given model and params and pool
   * \param m the model to operate on
   * \param p the params to use
   * \param tpool the already init'd thread pool
   */
  Hogwild(Model &m, Params &p, hazy::thread::ThreadPool &tpool) :
      model_(m), params_(p), tpool_(tpool) { }

  /*! \brief Invoke Exec::Execute on just them odel and params in each thread
   * Do not use any examples, just parallelize with each thread working
   * on the model and parameters
   */
  template <class Exec>
  void WildExecute();

  /*! \brief Split the contents of the scanner's examples between each thread
   * Invokes Exec::Execute(model, params, example) for each example
   * Uses the permutation given by the scanner to determine how examples
   * are split between threads.
   */
  template <class Exec, template <class E> class Scan, class Ex>
  void ForEach(Scan<Ex> &scan);

  /*! \brief Split the contents of the scanner's examples between each thread
   * Like ForEach, except this also aggregates the result and returns it.
   */
  template <class Exec, template <class E> class Scan, class Ex>
  typename Exec::Aggregate_t AggForEach(Scan<Ex> &scan);

  /*! \brief Split the contents of the scanner between each thread by the block
   * Each thread is given a block of examples and it is up to the individual
   * threads how the work is split up. Aggregates the result from each thread
   * and returns it.
   */
  template <class Exec, template <class E> class Scan, class Ex>
  typename Exec::Aggregate_t AggBlockScan(Scan<Ex> &scan);
 
 private:
  Model &model_; //!< the model
  Params &params_; //!< the params
  hazy::thread::ThreadPool &tpool_; //!< thread pool to use for train & test
};

} // namespace hogwild
} // namespace hazy
#endif
