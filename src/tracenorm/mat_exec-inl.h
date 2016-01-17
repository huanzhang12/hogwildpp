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

#ifndef HAZY_TRACENORM_MAT_EXEC_INL_H
#define HAZY_TRACENORM_MAT_EXEC_INL_H

#include "hazy/vector/dot-inl.h"
#include "hazy/vector/scale_add-inl.h"
#include "hazy/vector/operations-inl.h"

#include "hazy/hogwild/tools-inl.h"


#include "mat_exec.h"

namespace hazy {
namespace hogwild {
namespace tnorm {

double ComputeLoss(MFModel const &m, types::Entry const &e) {
  float err = vector::Dot(m.L[e.row],m.R[e.col]) + m.mean - e.rating;
  return err*err;
}

inline void ModelUpdate(MFModel &model, MFParams const &params, 
                        types::Entry const &e, 
                        vector::FVector<double> &swapL) {
  // declare things used later
  int row_index = e.row;
  int col_index = e.col;
  double mu = params.mu;
  double stepsize = params.step_size;
  int const * const L_degree = params.L_degree;
  int const * const R_degree = params.R_degree;

  vector::FVector<double> &Li = model.L[row_index];
  vector::FVector<double> &Rj = model.R[col_index];

  // Compute the error
  double err = vector::Dot(Li,Rj) + params.mean - e.rating;
  double er = -(params.step_size * err);

  double scal = 1 - params.mu*params.step_size/((double) L_degree[row_index]);

  // update the model
  vector::ScaleInto(Li, scal, swapL);
  vector::ScaleAndAdd(swapL, Rj, er);
  assert(R_degree[col_index] != 0);
  float scal_r = 1 - mu*stepsize/((float) R_degree[col_index]);
  vector::Scale(Rj, scal_r);
  vector::ScaleAndAdd(Rj, Li, er);

  // Store back L
  vector::CopyInto(swapL, Li);
}


double MFExec::UpdateModel(TNormTask &task, unsigned tid, unsigned total) {
  MFModel &model = *task.model;
  MFParams const &params = *task.params; 
  vector::FVector<types::Entry> const & exampsvec = task.block->ex;

  assert(params.L_degree != params.R_degree);
  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);

  types::Entry const * const examps = exampsvec.values;
  // declare some swap space on the stack for the update rule
  double swap[params.max_rank+1];
  vector::FVector<double> swapL(swap, params.nRows);

  size_t s = 0;
  // individually update the model for each example
  for (unsigned i = start; i < end; i++) {
    ModelUpdate(model, params, examps[i], swapL);
    s++;
  }
  return 0.0;
}

double MFExec::TestModel(TNormTask &task, unsigned tid, unsigned total) {
  MFModel const &model = *task.model;
  //MFParams const &params = *task.params; 
  vector::FVector<types::Entry> const & exampsvec = task.block->ex;

  // calculate which chunk of examples we work on
  size_t start = hogwild::GetStartIndex(exampsvec.size, tid, total); 
  size_t end = hogwild::GetEndIndex(exampsvec.size, tid, total);

  // keep const correctness
  types::Entry const * const examps = exampsvec.values;
  double loss = 0.0;
  // compute the loss for each example
  for (unsigned i = start; i < end; i++) {
    loss += ComputeLoss(model, examps[i]);
  }
  // return the number of examples we used and the sum of the loss
  return loss;
} 

} // namespace matfact
} // namespace hogwild

} // namespace hazy
#endif

