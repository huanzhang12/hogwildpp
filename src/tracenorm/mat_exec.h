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

#ifndef HAZY_HOGWILD_MATFACT_MAT_EXEC_H
#define HAZY_HOGWILD_MATFACT_MAT_EXEC_H
#include <iostream>
#include <cstdio>

#include "hazy/types/tuple.h"

namespace hazy {
namespace hogwild {
namespace tnorm {

//! Modifies the model using the single example
void inline ModelUpdate(MFModel &model, MFParams const &params, 
                        types::Entry const &e);
    
//! Returns the loss of the single example using the given model
double ComputeLoss(const MFModel &m, const types::Entry &e);

//! Container for methods to train and test an SVM
class MFExec {
 public:
  /// Preforms updates to the model
  /*! Updates by scanning over examples, uses the thread id and total
   * number of threads to determine which chunk of examples to work on.
   * \param task container of model, params, and examples 
   * \param tid the thread ID; 0 <= tid < total
   * \param total the total number of threads working on updating
   */
  static double UpdateModel(TNormTask &task, unsigned tid, unsigned total);

  /// Compute error of the task's model and the task's examples
  /*! Computes the error of each example with given the model. Uses the
   * number of threads to determine which chunk of examples to work on.
   * TODO XXX Needs to aggregate the RMSE
   * \param task container of model, params, and examples
   * \param tid the thread ID; 0 <= tid < total
   * \param total the total number of threads working on updating
   */
  static double TestModel(TNormTask &task, unsigned tid, unsigned total);

  //! Invoked after each training epoch, causes the stepsize to decay
  static void PostUpdate(MFModel &model, MFParams &params) {
    params.step_size *= params.step_decay;
  }
  static void PostEpoch(MFModel &model, MFParams &params) {
  }
};

} // namespace matfact
} // namespace hogwild
} // namespace hazy
#endif
