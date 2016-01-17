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

#include <cstdlib>

#include "hazy/hogwild/hogwild-inl.h"
#include "hazy/hogwild/memory_scan.h"
#include "hazy/hogwild/file_scan.h"
#include "hazy/scan/tsvfscan.h"
#include "hazy/scan/binfscan.h"

#include "mat_model.h"
#include "mat_exec-inl.h"

// Hazy imports
using namespace hazy;
using namespace hazy::hogwild;
using scan::TSVFileScanner;
using scan::MatlabTSVFileScanner;

using namespace hazy::hogwild::tnorm;

//using namespace hazy::hogwild::tnorm;




int main(int argc, char* argv[]) {

  if (argc != 8) {
    printf("usage: predict <model prefix> <mean> <nrows> <ncols> <rank> <entries to predict> <output file>\n");
    printf("model prefix: path to prefix-L.txv and prefix-R.tsv\n");
    printf("mean: float, the mean.\n");
    printf("nrows: number of rows in the matrix (rows in L)\n");
    printf("ncols: number of columns in the matrix (columns in R)\n");
    printf("rank: number of columns in L, rows in R\n");
    printf("entries to predict: path to a .tsv of pairs of integers indicating which entry to predict (base 0 index!)\n");
    printf("output file: the file to output\n");
    return -1;
  }

  const char *fileprefix = argv[1];
  double mean = atof(argv[2]);
  size_t rows = atol(argv[3]);
  size_t cols = atol(argv[4]);
  size_t rank = atol(argv[5]);
  const char *topredict_file = argv[6];
  const char *out_file = argv[7];

  // keep track of the number of non-zeros in each row
  int *row_counts = new int[rows];
  // keep track of the numbero f non-zeros in each column
  int *col_counts = new int[cols];

  printf("Loading from prefix `%s'\n", fileprefix);
  printf("#rows = %lu, #cols = %lu, rank = %lu\n", rows, cols, rank);
  printf("Output will have %lu entries!\n", rows * cols);

  Model_t model (mean, rows, cols, rank);

  char buf[1024];
  sprintf(buf, "%s-L.tsv", fileprefix);
  scan::TSVFileScanner lscan(buf);
  while (lscan.HasNext()) {
    const types::Entry &e = lscan.Next();
    row_counts[e.row]++;
    model.L[e.row].values[e.col] = e.rating;
  }

  sprintf(buf, "%s-R.tsv", fileprefix);
  printf("Loading from %s...\n", buf);
  scan::TSVFileScanner rscan(buf);
  while (rscan.HasNext()) {
    const types::Entry &e = rscan.Next();
    col_counts[e.row]++;
    model.R[e.row].values[e.col] = e.rating;
  }

  printf("Loaded the model!\n");

  FILE *out = fopen(out_file, "w");

  scan::OffsetTSVFileScanner<0, 1> predict_scan(topredict_file);
  size_t count = 0;
  while (predict_scan.HasNext()) {
    count++;
    const types::Tuple<1> &t = predict_scan.Next();
    int row = t.posn[0];
    int col = (int) t.rating;
    assert(row >= 0);
    assert(col >= 0);

    // if the row or column is not present, predict the mean
    if ((row_counts[row] == 0) || (col_counts[col] == 0)) {
      printf("Predicting mean!\n");
      fprintf(out, "%d\t%d\t%f\n", row, col, mean);
    } else {
      fprintf(out, "%d\t%d\t%f\n", row, col, vector::Dot(model.L[row], model.R[col])+mean);
    }
  }


  printf("Predicted %lu entries\n", count);
  return 0;
}
