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

#ifndef HAZY_HOGWILD_FREEFORALL_H
#define HAZY_HOGWILD_FREEFORALL_H

#include "hazy/vector/fvector.h"
#include "hazy/thread/thread_pool-inl.h"

#include "hazy/hogwild/hogwild_task.h"

namespace hazy {
namespace hogwild {

/*! Executes the hogwild task using in the thread pool, returns the results
 * Each thread in the pool will return a double, these are returned in
 * result which is to be already allocated to be the right size (nthreads)
 * \param task the howgwild task to execute
 * \param tpool the thread pool to use
 * \param hook the update function to run in each thread on the task
 * \param result the array to ADD the result of each thread to
 */
template <class HogwildTask_t>
void FreeForAll(HogwildTask_t &task, hazy::thread::ThreadPool &tpool,
                double (*hook)(HogwildTask_t&, unsigned, unsigned),
                vector::FVector<double> &result);


/*! 
 * Each thread in the pool will return a double, these are returned in
 * result which is to be already allocated to be the right size (nthreads)
 * \param m the model to use
 * \param p the params to use
 * \param scan the scanner which implements ExampleBlock<Example>& Next()
 * \param tpool the thread pool for execution
 * \param hook the update function to run in each thread on the task
 * \param result the array to ADD the result of each call to hook, based on
 *    on the thread id it was executed with.
 */
template <class Model, class Params, class Example, class Scan>
size_t FFAScan(Model &m, Params &p, Scan &scan, 
               hazy::thread::ThreadPool &tpool,
     double (*hook)(HogwildTask<Model, Params, Example>&, unsigned, unsigned),
               vector::FVector<double> &result);

template <class Model, class Params, class Example, class Scan>
size_t FFAScan(Model &m, Params &p, Scan &scan, 
               hazy::thread::ThreadPool &tpool,
     double (*hook)(HogwildTask<Model, Params, Example>&, unsigned, unsigned),
               vector::FVector<double> &result, 
               hazy::util::Clock &train_time, hazy::util::Clock &epoch_time); 
} // namespace hogwild
} // namespace hazy
#endif
