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

#include "hazy/vector/dot-inl.h"
#include "hazy/vector/scale_add-inl.h"
#include "hazy/hogwild/tools-inl.h"

namespace hazy {
namespace hogwild {
namespace svm {

fp_type inline ComputeLoss(const SVMExample &e, const SVMModel& model) {
  // determine how far off our model is for this example
  vector::FVector<fp_type> const &w = model.weights;
  fp_type dot = vector::Dot(w, e.vector);
  return std::max(1 - dot * e.value, static_cast<fp_type>(0.0));
}

int inline ComputeAccuracy(const SVMExample &e, const SVMModel& model) {
  // determine how far off our model is for this example
  vector::FVector<fp_type> const &w = model.weights;
  fp_type dot = vector::Dot(w, e.vector);
  return !!std::max(dot * e.value, static_cast<fp_type>(0.0));
}

void inline ModelUpdate(const SVMExample &examp, const SVMParams &params, 
                 SVMModel *model) {
  vector::FVector<fp_type> &w = model->weights;

  // evaluate this example
  fp_type wxy = vector::Dot(w, examp.vector);
  wxy = wxy * examp.value;

  if (wxy < 1) { // hinge is active.
    fp_type const e = params.step_size * examp.value;
    vector::ScaleAndAdd(w, examp.vector, e);
  }

  fp_type * const vals = w.values;
  unsigned const * const degs = params.degrees;
  size_t const size = examp.vector.size;

  // update based on the evaluation
  fp_type const scalar = params.step_size * params.mu;
  for (int i = size; i-- > 0; ) {
    int const j = examp.vector.index[i];
    unsigned const deg = degs[j];
    vals[j] *= 1 - scalar / deg;
  }
}

void SVMExec::PostUpdate(SVMModel &model, SVMParams &params) {
  // Reduce the step size to encourage convergence
  params.step_size *= params.step_decay;
}

double SVMExec::UpdateModel(SVMTask &task, unsigned tid, unsigned total) {

  SVMModel  &model = *task.model;

  SVMParams const &params = *task.params;
  vector::FVector<SVMExample> const & exampsvec = task.block->ex;
  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);
  // optimize for const pointers 
  size_t *perm = task.block->perm.values;
  SVMExample const * const examps = exampsvec.values;
  SVMModel * const m = &model;
  // individually update the model for each example
  // printf("UpdateModel: thread id %d updating model from %lu to %lu\n", tid, start, end);
  for (unsigned i = start; i < end; i++) {
    size_t indirect = perm[i];
    ModelUpdate(examps[indirect], params, m);
  }
  return 0.0;
}

double SVMExec::TestModel(SVMTask &task, unsigned tid, unsigned total) {
  SVMModel const &model = *task.model;

  //SVMParams const &params = *task.params;
  vector::FVector<SVMExample> const & exampsvec = task.block->ex;

  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);

  // keep const correctness
  SVMExample const * const examps = exampsvec.values;
  fp_type loss = 0.0;
  // compute the loss for each example
  for (unsigned i = start; i < end; i++) {
    fp_type l = ComputeLoss(examps[i], model);
    loss += l;
  }
  // return the number of examples we used and the sum of the loss
  //counted = end-start;
  return loss;
}

double SVMExec::ModelAccuracy(SVMTask &task, unsigned tid, unsigned total) {
  SVMModel const &model = *task.model;

  //SVMParams const &params = *task.params;
  vector::FVector<SVMExample> const & exampsvec = task.block->ex;

  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);

  // keep const correctness
  SVMExample const * const examps = exampsvec.values;
  int correct = 0;
  // compute the loss for each example
  for (unsigned i = start; i < end; i++) {
    int l = ComputeAccuracy(examps[i], model);
    correct += l;
  }
  // return the number of examples we used and the sum of the loss
  //counted = end-start;
  return correct;
}

double SVMExec::ModelObj(SVMTask &task, unsigned tid, unsigned total) {
  SVMModel const &model = *task.model;

  //SVMParams const &params = *task.params;
  vector::FVector<SVMExample> const & exampsvec = task.block->ex;

  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);

  // keep const correctness
  SVMExample const * const examps = exampsvec.values;
  fp_type loss = 0.0;
  // compute the loss for each example
  for (unsigned i = start; i < end; i++) {
    fp_type l = ComputeLoss(examps[i], model);
    loss += l;
  }
  start = hogwild::GetStartIndex(model.weights.size, tid, total);
  end = hogwild::GetEndIndex(model.weights.size, tid, total);
  double const * const weights = model.weights.values;
  fp_type reg = 0.0;
  // compute the regularization term
  for (unsigned i = start; i < end; ++i) {
    reg += weights[i] * weights[i];
  }
  // return the number of examples we used and the sum of the loss
  //counted = end-start;
  return loss + 0.5 * reg;
}

} // namespace svm
} // namespace hogwild

} // namespace hazy
