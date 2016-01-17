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


#ifndef HAZY_HOGWILD_INSTANCES_MATFACT_MAT_MODEL_H
#define HAZY_HOGWILD_INSTANCES_MATFACT_MAT_MODEL_H
#include <iostream>
#include <cstdio>
#include <cstdlib>

#include "hazy/vector/fvector.h"
#include "hazy/vector/svector.h"
#include "hazy/scan/tsvfscan.h"

#include "hazy/hogwild/hogwild_task.h"

namespace hazy {
namespace hogwild {
//! Matrix Factorization
namespace tnorm {

//! Parameters for Tracenorm
struct MFParams {
  int *L_degree;
  int *R_degree;
  size_t nExamples, max_rank;
  double mean, mu;

  double step_size;
  double step_decay;

  size_t nRows, nCols;

  //! Expects to have the map applied below, parameter_map
  void Setup(size_t _nRows, size_t _nCols, size_t _nExamples) { 
    assert(_nRows > 0 );
    assert(_nCols > 0 );
    nRows = _nRows;
    nCols = _nCols;
    nExamples = _nExamples;
    L_degree = new int[nRows]; memset(L_degree, 0, sizeof(int)*nRows);
    R_degree = new int[nCols]; memset(R_degree, 0, sizeof(int)*nCols);
    mean     = 0.0;
  }

  explicit MFParams(const MFParams &src) {
    L_degree = src.L_degree;
    R_degree = src.R_degree;
    nExamples = src.nExamples;
    max_rank = src.max_rank;
    mean = src.mean;
    mu = src.mu;
    step_size = src.step_size;
    step_decay = src.step_decay;
    nRows = src.nRows;
    nCols = src.nCols;
  }

  //! Defalut constructor
  explicit MFParams() : L_degree(NULL), R_degree(NULL), nExamples(-1) { }

};

//! Tracenorm Model
struct MFModel {
  vector::FVector<double> *L; //!< Left matrix
  vector::FVector<double> *R; //!< Right matrix
  double mean;
  size_t rows, cols;
  size_t max_rank;

  explicit MFModel() { }
  //! Allocates the left and right matrix
  explicit MFModel(double mean, size_t rows, size_t cols, size_t max_rank) {
    // Allocate and init the matrix
    L = new vector::FVector<double>[rows];
    this->rows = rows;
    this->cols = cols;
    this->max_rank = max_rank;
    this->mean = mean;
    for (size_t i = 0; i < rows; i++) { //p.nRows; i++) {
      double *arr = new double[max_rank];
      vector::FVector<double> &v = L[i];
      v.size = max_rank;
      v.values = arr;
      for(size_t j = 0; j < max_rank; j++) {
        // the magic initial values....
        L[i].values[j] = (drand48()-0.5)*1e-3;
      }
    }

    // Allocate the matrix
    R = new vector::FVector<double>[cols];
    for (size_t i = 0; i < cols; i++) {
      R[i].size = max_rank;
      R[i].values = new double[max_rank];
      for(size_t j = 0; j < max_rank; j++) {
        // chosen by magic
        R[i].values[j] = (drand48()-0.5)*1e-3;
      }
    }
    printf("Mean = %f\n", mean);
  }

  MFModel* Clone() {
    MFModel *p = new MFModel();
    p->rows = rows;
    p->cols = cols;
    p->mean = mean;
    p->max_rank = max_rank;
    p->L = new vector::FVector<double>[rows];
    for (size_t i = 0; i < rows; i++) { //p.nRows; i++) {
      vector::FVector<double> &v = p->L[i];
      v.values = new double[max_rank];
      v.size = L[0].size;
    }
    p->R = new vector::FVector<double>[cols];
    for (size_t i = 0; i < cols; i++) { //p.nRows; i++) {
      vector::FVector<double> &v = p->R[i];
      v.values = new double[max_rank];
      v.size = R[0].size;
    }
    p->CopyFrom(*this);
    return p;
  }

  void CopyFrom(MFModel &m) {
    for (size_t i = 0; i < rows; i++) { //p.nRows; i++) {
      memcpy(L[i].values, m.L[i].values, 
             sizeof(double)*L[0].size);
    }
    for (size_t i = 0; i < cols; i++) { //p.nRows; i++) {
      memcpy(R[i].values, m.R[i].values, 
             sizeof(double)*R[0].size);
    }
  }

  void OutputToFile(char const *name, const MFParams &p) {
    char buf[1024];
    sprintf(buf, "%s-L.tsv", name);
    printf("Writing to %s\n", buf);
    FILE* out = fopen(buf, "w");
    for (unsigned i = 0; i < p.nRows; i++) {
      for (unsigned j = 0; j < L[i].size; j++) {
        fprintf(out, "%u\t%u\t%f\n", i, j, L[i].values[j]);
      }
    }
    fclose(out);
    sprintf(buf, "%s-R.tsv", name);
    printf("Writing to %s\n", buf);
    out = fopen(buf, "w");
    for (unsigned i = 0; i < p.nCols; i++) {
      for (unsigned j = 0; j < R[i].size; j++) {
        fprintf(out, "%u\t%u\t%f\n", i, j, R[i].values[j]);
      }
    }
    fclose(out);
  }

  void LoadFromFile(char const *name) {
    char buf[1024];
    sprintf(buf, "%s-L.tsv", name);
    printf("Loading from %s...\n", buf);
    scan::TSVFileScanner lscan(buf);
    while (lscan.HasNext()) {
      const types::Entry &e = lscan.Next();
      L[e.row].values[e.col] = e.rating;
    }

    sprintf(buf, "%s-R.tsv", name);
    printf("Loading from %s...\n", buf);
    scan::TSVFileScanner rscan(buf);
    while (rscan.HasNext()) {
      const types::Entry &e = rscan.Next();
      R[e.row].values[e.col] = e.rating;
    }
  }
};


typedef MFParams Params_t;
typedef MFModel Model_t;


//! Adds the entrie's row and column to the degree counts, and updates mean
void parameter_map(MFParams &tnp, types::Entry const &e) {
  tnp.L_degree[e.row] ++;
  tnp.R_degree[e.col] ++;
  tnp.mean += e.rating/((double) tnp.nExamples);
}
typedef HogwildTask<Model_t, Params_t, types::Entry> TNormTask;

template <class Scan>
void SetParamsByScan(Scan &scan, Params_t &params) {
  size_t nex = 0;
  size_t max_row = 0;
  size_t max_col = 0;
  
  scan.Reset();
  while (scan.HasNext()) {
    const types::Entry &e = scan.Next();
    if ( (e.col > 0) && (static_cast<size_t>(e.col) > max_col) ) 
      max_col = e.col;
    if ( (e.row > 0) && (static_cast<size_t>(e.row) > max_row) ) 
      max_row = e.row;
    nex++;
  }
  params.Setup(max_row+1, max_col+1, nex);

  scan.Reset();
  while (scan.HasNext()) {
    const types::Entry &e = scan.Next();
    parameter_map(params, e);
  }
}

} // namespace matfact
} // namespace hogwild
} // namespace hazy

#endif
