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

#ifndef HAZY_HOGWILD_INSTANCES_CUTS_CUT_EXEC_H
#define HAZY_HOGWILD_INSTANCES_CUTS_CUT_EXEC_H
#include <set>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "hazy/vector/fvector.h"
#include "hazy/types/tuple.h"
#include "cut_model.h"

namespace hazy {
namespace hogwild {
namespace cuts {

std::set<int>
load_intlist(char *fname) {
  std::set<int> s;
  std::ifstream f_in;
  int i;
  f_in.open(fname, std::ios::in);
  while(!f_in.eof()) {
    f_in >> i;
    s.insert(i - 1);
  }
  f_in.close();
  return s;
}

//! Clips to between [0, 1]
double clip(double v) { return std::min(std::max(v, 0.0),1.0); }
//1 Returns the sign (as 1 or -1) of v
double sign(double v) { if(v == 0.0) return 0.0; 
  return (v > 0.0) ? 1.0 : -1.0; }

//! Updates the model using the given example
/*! This is a razy operation that will cause terrible cache behavior, most
 * likely.
 */
void ModelUpdate(types::Entry const &e, CutParams const &params, 
                 CutModel &model) {
  int row_index = e.row;
  int col_index = e.col;
  float weight = e.rating;
  vector::FVector<double> &row_model = model.get(row_index);
  vector::FVector<double> &col_model = model.get(col_index);
  for(int d = model.get_dim() - 1; d >= 0; d--) {
    double   x = row_model.values[d];
    double   y = col_model.values[d];
    double err = weight*sign(x - y);
    row_model.values[d] = x - params.stepsize*err;
    col_model.values[d] = y + params.stepsize*err;     
  }
  vector::SimplexProject(row_model);
  vector::SimplexProject(col_model);
}

//! Computes the loss for the single model & example pair
double ComputeLoss(CutModel const &m, types::Entry const &e) {
  int dmax =  m.get_dim();
  int row_index = e.row;
  int col_index = e.col;
  double weight = e.rating;

  double err_total = 0.0;
  for (int d = 0; d < dmax; d++) {
    double x = m.get(row_index, d);
    double y = m.get(col_index, d);
    err_total += weight*std::abs(x - y);
  }
  return err_total;
}

//! Computes the loss for the single model & example pair, as 0 or 1
float ZeroOneLoss(CutModel const &m, types::Entry const &e) {
  int dmax         = m.get_dim();
  int row_index    = e.row;
  int col_index    = e.col;
  double weight    = e.rating;

  int row_max_index = 0, col_max_index  = 0;
  double row_max = m.get(row_index,0);
  double col_max = m.get(col_index,0);
  for(int d = 0; d < dmax; d++) {
    double rv = m.get(row_index, d);
    double cv = m.get(col_index, d);
    if(rv > row_max) { 
      row_max_index = d ;
      row_max       = rv;
    }
    if(cv > col_max) {
      col_max_index = d;
      col_max       = cv;
    }
  }

  double err_total = 0.0;
  if(row_max_index != col_max_index) {
    err_total = weight;
  }
  return err_total;
}

//! Container for executing Hogwild! operations
/*! Exports the 3 functions needed by Hogwild! FreeForAll to train the model.
 * \implements Exec
 */
class CutExec {
 public:
  //! Flag set if ZeroOneLoss is to be used; otherwise uses regular loss
  static bool UseZeroOneLoss;
  static double UpdateModel(CutTask &task,  unsigned tid, unsigned total) {
    CutModel &model = *task.model;
    CutParams const &params = *task.params;
    vector::FVector<types::Entry> const &exampsvec = task.block->ex;
    // calculate which chunk of examples we work on
    unsigned start = (exampsvec.size / total) * tid;
    unsigned end = std::min((exampsvec.size / total) * (tid + 1),
                       exampsvec.size);
    if (tid == total-1) {
      end = exampsvec.size;
    }
    types::Entry const * const examps = exampsvec.values;
    CutModel * const m = &model;
    // individually update the model for each example
    for (unsigned i = start; i < end; i++) {
      ModelUpdate(examps[i], params, *m);
    }
    return 0.0;
  }

  static double TestModel(CutTask &task, unsigned tid, unsigned total) {
    CutModel const &model = *task.model;
    vector::FVector<types::Entry> const &exampsvec = task.block->ex;
    // calculate which chunk of examples we work on
    unsigned start = (exampsvec.size / total) * tid;
    unsigned end = std::min((exampsvec.size / total) * (tid + 1),
                       exampsvec.size);
    if (tid == total-1) {
      end = exampsvec.size;
    }

    types::Entry const * const examps = exampsvec.values;
    double loss = 0.0;
    // compute the loss for each example
    for (unsigned i = start; i < end; i++) {
      double l;
      if (UseZeroOneLoss) {
        l = ZeroOneLoss(model, examps[i]);
      } else {
        l = ComputeLoss(model, examps[i]);
      }
      loss += l;
    }
    return loss;
  }

  static void PostUpdate(CutModel &m, CutParams &p) {
    p.stepsize *= p.step_diminish;
  }
  static void PostEpoch(CutModel &m, CutParams &p) {
    p.stepsize *= p.step_diminish;
  }
};

bool CutExec::UseZeroOneLoss = false;

} // namesace cuts
} // namespace hogwild

} //namespace hazy
#endif
