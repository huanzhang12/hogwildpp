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

#ifndef HAZY_HOGWILD_INSTANCES_SVM_SVM_EXEC_H
#define HAZY_HOGWILD_INSTANCES_SVM_SVM_EXEC_H

#include <cmath>

#include "hazy/hogwild/hogwild_task.h"

#include "svmmodel.h"

namespace hazy {
namespace hogwild {
namespace svm {

//! Changes the model using the given example 
void inline ModelUpdate(const SVMExample &examp, const SVMParams &params, 
                 NumaSVMModel *model, size_t &updates, size_t &count);

//! Returns the loss for the given example and model
fp_type inline ComputeLoss(const SVMExample &e, const NumaSVMModel& model);

//! Container for methods to train and test an SVM
class NumaSVMExec {
 public:
  /// Preforms updates to the model
  /*! Updates by scanning over examples, uses the thread id and total
   * number of threads to determine which chunk of examples to work on.
   * \param task container of model, params, and examples 
   * \param tid the thread ID; 0 <= tid < total
   * \param total the total number of threads working on updating
   */
  static double UpdateModel(SVMTask &task, unsigned tid, unsigned total);

  /// Compute error of the task's model and the task's examples
  /*! Computes the error of each example with given the model. Uses the
   * number of threads to determine which chunk of examples to work on.
   * TODO XXX Needs to aggregate the RMSE
   * \param task container of model, params, and examples
   * \param tid the thread ID; 0 <= tid < total
   * \param total the total number of threads working on updating
   */
  static double TestModel(SVMTask &task, unsigned tid, unsigned total);

  //! Invoked after each training epoch, causes the stepsize to decay
  static void PostUpdate(NumaSVMModel &model, SVMParams &params);

  static void PostEpoch(NumaSVMModel &model, SVMParams &params) {
  }
  static double ModelObj(SVMTask &task, unsigned tid, unsigned total);
  static double ModelAccuracy(SVMTask &task, unsigned tid, unsigned total);
 private:
  static int GetNumaNode();
  static int GetLatestModel(SVMTask &task, unsigned tid, unsigned total);
};

} // namespace svm
} // namespace hogwild

} // namespace hazy
#include "svm_exec.hxx"
#endif
