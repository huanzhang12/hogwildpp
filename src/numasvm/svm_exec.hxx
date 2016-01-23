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
                 NumaSVMModel *model, NumaSVMModel *next_model, int tid) {
  vector::FVector<fp_type> &w = model->weights;
  // vector::FVector<fp_type> &dw = model->delta_weights;

  // evaluate this example
  fp_type wxy = vector::Dot(w, examp.vector);
  // fp_type wxy = vector::AddAndDot(w, dw, examp.vector);
  wxy = wxy * examp.value;

  if (wxy < 1) { // hinge is active.
    fp_type const e = params.step_size * examp.value;
    vector::ScaleAndAdd(w, examp.vector, e);
    // vector::ScaleAndAdd(dw, examp.vector, e);
  }

  fp_type * const vals = w.values;
  // fp_type * const dvals = dw.values;
  unsigned const * const degs = params.degrees;
  size_t const size = examp.vector.size;
  // update based on the evaluation
  fp_type const scalar = params.step_size * params.mu;
  for (int i = size; i-- > 0; ) {
    int const j = examp.vector.index[i];
    unsigned const deg = degs[j];
    vals[j] *= (1 - scalar / deg);
    // dvals[j] = dvals[j] * (1 - scalar / deg) - vals[j] * (scalar / deg);
  }
  // update dw to the other thread
  if (model->GetAtomic() == tid) { // TODO: this should be physical cpu number
    if (next_model) {
      fp_type * const old_vals = model->delta_weights.values;
      fp_type * const next_vals = next_model->weights.values;
      fp_type * const next_old_vals = next_model->delta_weights.values;
      int counter = 0;
      for (unsigned i = 0; i < w.size; ++i) {
        fp_type wi = vals[i];
        fp_type delta = wi - old_vals[i];
        fp_type next = next_vals[i];
        if (fabs(delta) > 1e-1) {
	  fp_type new_wi = next + delta;
	  next_vals[i] = new_wi;
	  vals[i] = new_wi;
          old_vals[i] = new_wi;
          next_old_vals[i] += delta;
	  // dvals[i] = 0;
          counter++;
        }
        else {
	  vals[i] = next + delta;
          old_vals[i] = next;
        }
      }
      // printf("%d:%d/%ld\n", tid, counter, w.size);
      model->IncAtomic();
    }
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
  // individually update the model for each example
  int weights_index = model.thread_to_weights_mapping[tid];
  int next_weights = model.next_weights[tid];
  NumaSVMModel * const m = &task.model[weights_index];
  NumaSVMModel * const next_m = next_weights >= 0 ? &task.model[next_weights] : NULL;
  int atomic_inc_value = m->atomic_inc_value;
  int atomic_mask = m->atomic_mask;
  printf("UpdateModel: thread id %d on node %d using data %p from %lu to %lu,"
         "and model %d->%d at %p->%p, (atomic+%d) & %x\n", 
         tid, node, exampsvec.values, start, end, weights_index, next_weights,
         m->weights.values, next_weights >= 0 ? next_m->weights.values: NULL, 
         atomic_inc_value, atomic_mask);
  for (unsigned i = start; i < end; i++) {
    size_t indirect = perm[i];
    ModelUpdate(examps[indirect], params, m, next_m, tid);
  }
  
  // vector::ScaleAndAdd(m->weights, m->delta_weights, 1);
  
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
