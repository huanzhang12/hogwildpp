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

#ifndef HAZY_HOGWILD_INSTANCES_CUTS_CUT_MODEL_H
#define HAZY_HOGWILD_INSTANCES_CUTS_CUT_MODEL_H
#include <set>
#include <algorithm>

#include "hazy/types/tuple.h"
#include "hazy/vector/fvector.h"
#include "hazy/vector/operations-inl.h"
#include "hazy/hogwild/hogwild_task.h"

namespace hazy {
namespace hogwild {
namespace cuts {

//! Parameters for multicut
struct CutParams {
  double stepsize;  
  double step_diminish;
  CutParams(float ss, float sd) : stepsize(ss), step_diminish(sd) { }
}; 

//! mutable model for multicut
class CutModel {
 public:
  int n, dim;
  std::set<int> &Terminals;
  vector::FVector<double> *weights;

  void CopyFrom(const CutModel &m) {
    vector::CopyInto(*m.weights, *weights);
  }

  CutModel* Clone() {
    CutModel *m = new CutModel(n, Terminals);
    m->CopyFrom(*this);
    return m;
  }

  void SumWeights() {
    double t = 0.0;
    size_t c = 0;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < dim; j++) {
        t += weights[i].values[j] * (j+1) / (i+1);
        c++;
      }
    }
    printf("Counted: %llu\n", static_cast<long long unsigned>(c));
    printf("Sum of weights: %f\n", t);
  }

  CutModel(int _n, std::set<int> &_terminals)
      : n(_n), Terminals(_terminals) {
    dim = (int) Terminals.size();
    weights = new vector::FVector<double>[n];
    for (int i = 0; i < n; i++) {
      weights[i].size = dim;
      weights[i].values = new double[dim];
    }

    std::cout << "Beginning simplex projects" << std::endl;
    for(int j = 0; j < n; j++) {
      //printf("... attempting = %d\n", j);
      vector::SimplexProject(weights[j]);
      //printf("j = %d / %d\n", j, n);
    }

    std::cout << "Zeroing " << dim << " terminals" << std::endl;      
    int counter = 0;
    for(std::set<int>::iterator i = Terminals.begin(); i != Terminals.end(); i++) {
      vector::Zero(weights[*i]);
      weights[*i].values[counter] = 1.0;
      counter++;
    }
  }

  ~CutModel() {
    delete weights;
  }

  // Accessors
  void set(int i, int j, double w) {
    if(Terminals.find(i) == Terminals.end()) { weights[i].values[j] = w; }	
  }

  vector::FVector<double>& get(int i) { return weights[i]; }
  void set(int i, const vector::FVector<double> &x) { 
    weights[i].size = x.size;
    weights[i].values = x.values;
  }
  float get(int i, int j) const { 
    assert(i >= 0 && i < n); 
    //printf("&get = %llux\n", reinterpret_cast<size_t>(&weights[i].values[j]));
    return weights[i].values[j];
  }
  void project(int i) {
    if(Terminals.find(i) == Terminals.end())
      vector::SimplexProject(weights[i]); 
  }
  int get_dim() const { return dim; }  

};

double compute_loss(const CutModel &m, const types::Entry &e);  

typedef HogwildTask<CutModel, CutParams, types::Entry> CutTask;

} // namesace cuts
} // namespace hogwild

} //namespace hazy
#endif
