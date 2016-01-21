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

#include <numa.h>
#include <sched.h>
#include <cstdio>

namespace hazy {
namespace hogwild {
namespace svm {

fp_type inline ComputeLoss(const SVMExample &e, const NumaSVMModel& model) {
  // determine how far off our model is for this example
  vector::FVector<fp_type> const &w = model.weights;
  fp_type dot = vector::Dot(w, e.vector);
  return std::max(1 - dot * e.value, static_cast<fp_type>(0.0));
}

void inline ModelUpdate(const SVMExample &examp, const SVMParams &params, 
                 NumaSVMModel *model) {
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

void NumaSVMExec::PostUpdate(NumaSVMModel &model, SVMParams &params) {
  // Reduce the step size to encourage convergence
  params.step_size *= params.step_decay;
}

int NumaSVMExec::GetNumaNode() {
  int cpu = sched_getcpu();
  return numa_node_of_cpu(cpu);
}

double NumaSVMExec::UpdateModel(SVMTask &task, unsigned tid, unsigned total) {

  int node = GetNumaNode();
  // TODO: per core model vector 
  NumaSVMModel  &model = *task.model;

  SVMParams const &params = *task.params;
  // Select the example vector array based on current node
  vector::FVector<SVMExample> const & exampsvec = task.block[node].ex;
  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);
  // optimize for const pointers 
  // Seclect the pointers based on current node
  size_t *perm = task.block[node].perm.values;
  SVMExample const * const examps = exampsvec.values;
  NumaSVMModel * const m = &model;
  // individually update the model for each example
  printf("UpdateModel: thread id %d on node %d updating model %p, data %p from %lu to %lu\n", tid, node, m->weights.values, exampsvec.values, start, end);
  for (unsigned i = start; i < end; i++) {
    size_t indirect = perm[i];
    ModelUpdate(examps[indirect], params, m);
  }
  return 0.0;
}

double NumaSVMExec::TestModel(SVMTask &task, unsigned tid, unsigned total) {
  int node = GetNumaNode();
  NumaSVMModel const &model = *task.model;

  //SVMParams const &params = *task.params;
  // Select the example vector array based on current node
  vector::FVector<SVMExample> const & exampsvec = task.block[node].ex;

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

} // namespace svm
} // namespace hogwild

} // namespace hazy
