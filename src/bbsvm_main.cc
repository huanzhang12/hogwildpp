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

#include "hazy/scan/tsvfscan.h"
#include "hazy/scan/binfscan.h"
#include "hazy/hogwild/hogwild-inl.h"
#include "hazy/hogwild/memory_scan.h"
#include "hazy/hogwild/bestball-inl.h"

#include "frontend_util.h"
#include "loader-inl.h"

#include "svm/svmmodel.h"
#include "svm/svm_loader.h"
#include "svm/svm_exec.h"


// Hazy imports
using namespace hazy;
using namespace hazy::hogwild;
using scan::TSVFileScanner;
using scan::MatlabTSVFileScanner;

using hazy::hogwild::svm::fp_type;


using namespace hazy::hogwild::svm;


typedef BestBall<SVMModel, SVMParams, SVMExample, SVMExec> AutoTuner;

int main(int argc, char** argv) {
  hazy::util::Clock wall_clock;
  wall_clock.Start();
  //Benchmark::StartExperiment(argc, argv);

  bool matlab_tsv = false;
  char *step_file = NULL;
  char *mufile = NULL;
  bool loadBinary = false;
  unsigned nepochs = 20;
  unsigned nthreads = 1;
  float step_decay = 0.8;
  static struct extended_option long_options[] = {
    {"mu", required_argument, NULL, 'u', "the maxnorm"},
    {"epochs"    ,required_argument, NULL, 'e', "number of epochs (default is 20)"},
    {"stepinitial",required_argument, NULL, 'i', "intial stepsize (default is 5e-2)"},
    {"step_decay",required_argument, NULL, 'd', "stepsize decay per epoch (default is 0.8)"},
    {"seed", required_argument, NULL, 's', "random seed (o.w. selected by time, 0 is reserved)"},
    {"splits", required_argument, NULL, 'r', "number of threads (default is 1)"},
    //{"shufflers", required_argument, NULL, 'q', "number of shufflers"},
    {"binary", required_argument,NULL, 'v', "load the file in a binary fashion"},
    {"matlab-tsv", required_argument,NULL, 'm', "load TSVs indexing from 1 instead of 0"},
    {NULL,0,NULL,0,0} 
  };

  char usage_str[] = "<train file> <test file>";
  int c = 0, option_index = 0;
  option* opt_struct = convert_extended_options(long_options);
  while( (c = getopt_long(argc, argv, "", opt_struct, &option_index)) != -1) 
  {
    switch (c) { 
      case 'v':
        loadBinary = (atoi(optarg) != 0);
        break;
      case 'm':
        matlab_tsv = (atoi(optarg) != 0);
        break;
      case 'u':
        mufile = optarg;
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
      case 'r':
        nthreads = atoi(optarg);
        break;
      case ':':
      case '?':
        print_usage(long_options, argv[0], usage_str);
        exit(-1);
        break;
    }
  }

  char * szTestFile, *szExampleFile;
  
  if(optind == argc - 2) {
    szExampleFile = argv[optind];
    szTestFile  = argv[optind+1];
  } else {
    print_usage(long_options, argv[0], usage_str);
    exit(-1);
  }
  //fp_type buf[50];

  vector::FVector<SVMExample> train_examps;
  vector::FVector<SVMExample> test_examps;

  size_t nfeats;

  if (loadBinary) {
    printf("Loading binary file...\n");
    scan::BinaryFileScanner scan(szExampleFile);
    nfeats = LoadSVMExamples(scan, train_examps);
    printf("Loaded binary file!\n");
  } else if (matlab_tsv) {
    MatlabTSVFileScanner scan(szExampleFile);
    nfeats = LoadSVMExamples(scan, train_examps);
  } else {
    TSVFileScanner scan(szExampleFile);
    nfeats = LoadSVMExamples(scan, train_examps);
  }
  if (matlab_tsv) {
    MatlabTSVFileScanner scantest(szTestFile);
    LoadSVMExamples(scantest, test_examps);
  } else {
    TSVFileScanner scantest(szTestFile);
    LoadSVMExamples(scantest, test_examps);
  }

  unsigned *degs = new unsigned[nfeats];
  printf("Loaded %lu examples\n", nfeats);
  for (size_t i = 0; i < nfeats; i++) {
    degs[i] = 0;
  }
  CountDegrees(train_examps, degs);

  printf("Loaded %lu test examples.\n", test_examps.size);
  printf("Loaded %lu train examples.\n", train_examps.size);
  printf("number of features: %lu\n", nfeats);


  // Create the list of different parameters to use, i.e. the params for each
  // ball in the 'best ball' algorithm
  std::vector<double> steps = load_floatlist(step_file);
  std::vector<double> mus = load_floatlist(mufile);
  assert(steps.size() > 0);
  assert(steps.size() == mus.size());
  vector::FVector<SVMParams*> pars (new SVMParams*[steps.size()], steps.size());
  for (unsigned i = 0; i < steps.size(); i++ ){
    printf("Ball #%i = (step=%lf, mu=%lf)\n", i, steps[i], mus[i]);
    pars.values[i] = new SVMParams(steps[i], step_decay, mus[i]);
    pars.values[i]->degrees = degs;
    pars.values[i]->ndim = nfeats;
  }

  // Initial copy of the model, the BestBall AutoTuner will .Clone() this
  SVMModel m(nfeats);

  MemoryScan<SVMExample> mscan(train_examps);
  MemoryScan<SVMExample> tscan(test_examps);

  AutoTuner bb (m, pars, nthreads);

  bb.RunExperiment(nepochs, wall_clock, mscan, tscan);

  return 0;
}

