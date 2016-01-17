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
#include <iostream>
#include <fstream>
#include <set>

#include "hazy/scan/binfscan.h"
#include "hazy/scan/tsvfscan.h"
#include "hazy/util/clock.h"
#include "hazy/thread/thread_pool-inl.h"

#include "hazy/hogwild/memory_scan.h"
#include "hazy/hogwild/bestball-inl.h"

#include "cuts/cut_model.h"
#include "cuts/cut_exec.h"

#include "frontend_util.h"
#include "loader-inl.h"

// Hazy imports
using namespace hazy;
using namespace hazy::hogwild;
using namespace hazy::hogwild::cuts;
// Hogwild! imports

using namespace std;

typedef BestBall<CutModel, CutParams, types::Entry, CutExec> AutoTuner;

int main(int argc, char** argv) {
  hazy::util::Clock wall_clock;
  wall_clock.Start();

  char* step_file = NULL;
  bool matlab_tsv = false;
  bool loadBinary = false;
  unsigned nepochs = 20;
  bool zero_one_loss  = false;
  unsigned nthreads = 1;
  float step_decay = 0.8;
  static struct extended_option long_options[] = {
    {"epochs"    ,required_argument, NULL, 'e', "number of epochs (default is 20)"},
    {"stepinitial",required_argument, NULL, 'i', "intial stepsize (default is 5e-2)"},
    {"step_decay",required_argument, NULL, 'd', "stepsize decay per epoch (default is 0.8)"},
    {"seed", required_argument, NULL, 's', "random seed (o.w. selected by time, 0 is reserved)"},
    {"splits", required_argument, NULL, 'r', "number of threads (default is 1)"},
    {"zero_one", required_argument, NULL, 'z', "use zero one loss (rounded solution)"},
    //{"shufflers", required_argument, NULL, 'q', "number of shufflers"},
    {"binary", required_argument,NULL, 'v', "load the file in a binary fashion"},
    {"matlab-tsv", required_argument,NULL, 't', "load TSVs indexing from 1 instead of 0"},
    {NULL,0,NULL,0,0} 
  };
  char usage_str[] = "<train file> <test file>";
  int c = 0, option_index = 0;
  option* opt_struct = convert_extended_options(long_options);

  while( (c = getopt_long(argc, argv, "", opt_struct, &option_index)) != -1) {
    switch (c) {
      case 'v':
        loadBinary = (atoi(optarg) != 0);
        break;
      case 'e':
        nepochs = atoi(optarg);
        break;
      case 'i':
        step_file = optarg;
        break;
      case 'd':
        step_decay = atof(optarg);
        break;
      case 't':
        matlab_tsv = (atoi(optarg) != 0);
        break;
      case 'r':
        nthreads = atoi(optarg);
        break;
      case 'z':
        zero_one_loss = (atoi(optarg) != 0);
        break;
      case ':':
      case '?':
        print_usage(long_options, argv[0], usage_str);
        exit(-1);
        break;
    }
  }

  char *szExampleFile, *terminalFile;
  
  if(optind == argc - 2) {
    szExampleFile = argv[optind];
    terminalFile    = argv[optind+1];
  } else {
    print_usage(long_options, argv[0], usage_str);
    exit(-1);
  }

  CutExec::UseZeroOneLoss = zero_one_loss;

  std::set<int> terminals = cuts::load_intlist(terminalFile);
  std::vector<double> step_sizes = load_floatlist(step_file);

  hazy::vector::FVector<CutParams*> pars(new CutParams*[step_sizes.size()],
                                   step_sizes.size());
  for (unsigned i = 0; i < step_sizes.size(); i++) {
    printf("Ball #%i = (step=%lf)\n", i, step_sizes[i]);
    pars.values[i] = new CutParams(step_sizes[i], step_decay);
  }


  hazy::vector::FVector<types::Entry> train_examps;
  size_t ncols;

  if (loadBinary) {
    scan::BinaryFileScanner scan(szExampleFile);
    ncols = LoadExamples(scan, train_examps);
  } else if (matlab_tsv) {
    scan::MatlabTSVFileScanner scan(szExampleFile);
    ncols = LoadExamples(scan, train_examps);
  } else {
    scan::TSVFileScanner scan(szExampleFile);
    ncols = LoadExamples(scan, train_examps);
  }

  printf("Loaded %lu train examples.\n", train_examps.size);
  printf("number of cols: %lu\n", ncols);
  fflush(stdout);


  CutModel m(ncols+1, terminals);

  MemoryScan<types::Entry> mscan(train_examps);
  AutoTuner bb(m, pars, nthreads);
  bb.RunExperiment(nepochs, wall_clock, mscan);

  return 0;
}

