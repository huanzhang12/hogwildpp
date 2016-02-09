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


#ifndef HAZY_HOGWILD_HOGWILD_H
#define HAZY_HOGWILD_HOGWILD_H

#include "hazy/vector/fvector.h"
#include "hazy/thread/thread_pool-inl.h"

namespace hazy {
namespace hogwild {

/*! \brief Hogwild! parallel executor
 * \tparam Exec implements Exec::UpateModel, Exec::TestModel, Exec::PostUpdate
 */
template <class Model, class Params, class Exec>
class Hogwild {
 public:
  /*! \brief Sets up for the given model and params and pool
   * \param m the model to operate on
   * \param p the params to use
   * \param tpool the already init'd thread pool
   */
  Hogwild(Model &m, Params &p, hazy::thread::ThreadPool &tpool) :
      model_(m), params_(p), tpool_(tpool) {
    res_.size = tpool_.ThreadCount();
    res_.values = new double[res_.size];
  }

  ~Hogwild() {
    delete [] res_.values;
  }

  /*! \brief Uses the scanner to upate the model
   * Resets the scanner then iterates until the end, updating the model
   */
  template <class Scan>
  void UpdateModel(Scan &scan);

  /*! \brief returns the RMSE of the model using examples form the scanner.
   * \param scan the scanner to use for computing the RMSE
   */
  template <class Scan>
  double ComputeRMSE(Scan &scan);

 
  template <class Scan>
  double ComputeObj(Scan &scan); 

  template <class Scan>
  double ComputeAccuracy(Scan &scan); 

  /*! \brief Runs an experiment printing statistics
   * \param nepochs number of epochs to run for
   * \param wall_clock a clock that was already started
   * \param trscan the train example scanner
   * \param tescan the test eample scanner
   */
  template <class TrainScan, class TestScan>
  void RunExperiment(int nepochs, hazy::util::Clock &wall_clock, 
                     TrainScan &trscan, TestScan &tescan);
 
  /*! \brief Runs an experiment printing statistics
   * \param nepochs number of epochs to run for
   * \param wall_clock a clock that was already started
   * \param trscan the train example scanner
   */
  template <class TrainScan>
  void RunExperiment(int nepochs, hazy::util::Clock &wall_clock, 
                     TrainScan &trscan);
 
 private:
  hazy::util::Clock train_time_; //!< measures total time spent in UpdateModel
  hazy::util::Clock test_time_; //!< measures total time in spent in ComputeRMSE
  hazy::util::Clock epoch_time_; //!< measures the time spent in most recent Update
  Model &model_; //!< the model
  Params &params_; //!< the params
  hazy::thread::ThreadPool &tpool_; //!< thread pool to use for train & test
  vector::FVector<double> res_; //!< results (for computing RMSE) size=nthreads

  /*! set the res_ to be all zeros
   */
  void Zero() { for (unsigned i = 0; i < res_.size; i++) res_.values[i] = 0; }
};

} // namespace hogwild
} // namespace hazy

#endif
