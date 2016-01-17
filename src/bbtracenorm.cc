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

#include "frontend_util.h"
#include "loader-inl.h"

#include "tracenorm/mat_model.h"
#include "tracenorm/mat_exec-inl.h"

#include "hazy/hogwild/bestball-inl.h"
// Hazy imports
using namespace hazy;
using namespace hazy::hogwild;
using scan::TSVFileScanner;
using scan::MatlabTSVFileScanner;

using namespace hazy::hogwild::tnorm;

//using namespace hazy::hogwild::tnorm;

typedef BestBall<Model_t, Params_t, types::Entry, MFExec> AutoTuner;

int main(int argc, char* argv[]) {
  hazy::util::Clock wall_clock;
  wall_clock.Start();
  // defaults
  Params_t p;
  int nEpochs = 20; 
  int seed = 0;
  int nSplits = 1;

  char *mu_file = NULL;
  char *step_file = NULL;

  p.step_size = 1e-1;
  p.step_decay = 1.0;
  p.mu       = -1.0;
  p.max_rank = 30;

  bool matlab_tsv = false;

  bool loadBinary = false;
  size_t bFileScan = 0;
  char *testFile = NULL, *outputTestFile = NULL;
  char *infile = NULL;
  static struct extended_option long_options[] = {
    {"mu", required_argument, NULL, 'u', "the maxnorm"},
    {"maxrank", required_argument, NULL, 'm', "size of factorization of L and R"},
    {"epochs"    ,required_argument, NULL, 'e', "number of epochs (default is 20)"},
    {"stepinitial",required_argument, NULL, 'i', "intial stepsize (default is 5e-2)"},
    {"step_decay",required_argument, NULL, 'd', "stepsize decay per epoch (default is 0.8)"},
    {"seed", required_argument, NULL, 's', "random seed (o.w. selected by time, 0 is reserved)"},
    {"splits", required_argument, NULL, 'r', "number of threads"},
    {"binary", required_argument,NULL, 'v', "load the file in a binary fashion"},
    {"file_scan", required_argument,NULL, 'f', "load the file in a scan (binary)"},
    {"outfile", required_argument,NULL, 'o', "write out to NAME-L.tsv and NAME-R.tsv"},
    {"infile", required_argument,NULL, 'l', "load from the given NAME-L.tsv and NAME-R.tsv"},
    {"matlab-tsv", required_argument,NULL, 't', "load TSVs indexing from 1 instead of 0"},
    {NULL,0,NULL,0,0} 
  };

  char usage_str[] = "<train file> <test file>";
  int c = 0, option_index = 0;
  option* opt_struct = convert_extended_options(long_options);
  while( (c = getopt_long(argc, argv, "", opt_struct, &option_index)) != -1) 
    {
      switch (c) {      
      case 'f':
    bFileScan     = (size_t) (((double)atof(optarg)) 
                              * 1024ULL * 1024ULL * 1024ULL);
  break;
      case 'w':
  testFile = optarg;
  break;
      case 'o':
  outputTestFile = optarg;
  break;
      case 'l':
  infile = optarg;
  break;
      case 'v':
  loadBinary = (atoi(optarg) != 0);
  break;
      case 't':
  matlab_tsv = (atoi(optarg) != 0);
  break;
      case 'q':
  assert(false);
  break;
      case 'u':
  mu_file = optarg;
  break;
      case 'm':
  p.max_rank  = atoi(optarg);
  break;
      case 'e':
  nEpochs = atoi(optarg);
  break;
      case 'i':
  step_file = optarg;
  break;
      case 'd':
  p.step_decay  = atof(optarg);
  break;
      case 's':
  seed              = atol(optarg);
  break;
      case 'r':
  nSplits         = atoi(optarg);
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

  if(mu_file == NULL) {
    std::cout << "mu is required parameter for the tracenorm" << std::endl;
    print_usage(long_options, argv[0], usage_str);
    exit(-1);
  }

  if(step_file == NULL) {
    std::cout << "step_size is required parameter for the tracenorm" << std::endl;
    print_usage(long_options, argv[0], usage_str);
    exit(-1);
  }

  vector::FVector<types::Entry> train_examps;
  vector::FVector<types::Entry> test_examps;
  Params_t &params = p;

  if (bFileScan) {
    printf("Loading binary file...\n");
    scan::BinaryFileScanner scan(szExampleFile);
    SetParamsByScan(scan, params);
    printf("done loading params.\n");
  } else if (loadBinary) {
    printf("Loading binary file...\n");
    scan::BinaryFileScanner scan(szExampleFile);
    SetParamsByScan(scan, params);
    LoadExamples(scan, train_examps);
    printf("Loaded binary file!\n");
  } else if (matlab_tsv) {
    MatlabTSVFileScanner scan(szExampleFile);
    SetParamsByScan(scan, params);
    scan.Reset();
    LoadExamples(scan, train_examps);
  } else {
    TSVFileScanner scan(szExampleFile);
    SetParamsByScan(scan, params);
    scan.Reset();
    LoadExamples(scan, train_examps);
  }
  if (matlab_tsv) {
    MatlabTSVFileScanner scan(szTestFile);
    scan.Reset();
    LoadExamples(scan, test_examps);
  } else {
    TSVFileScanner scan(szTestFile);
    LoadExamples(scan, test_examps);
  }

  printf("Loaded %lu test examples.\n", test_examps.size);
  printf("Loaded %lu train examples.\n", params.nExamples);
  printf("number of rows: %lu\n", params.nRows);
  printf("number of cols: %lu\n", params.nCols);
  fflush(stdout);

  std::vector<double> step_sizes = load_floatlist(step_file);
  std::vector<double> mus = load_floatlist(mu_file);

  assert(step_sizes.size() > 0);
  assert(step_sizes.size() == mus.size());
  for (unsigned i = 0; i < mus.size(); i++ ){
    printf("Ball #%i = (step=%lf, mu=%lf)\n", i, step_sizes[i], mus[i]);
  }

  //int nmodels = 2;
  Model_t model (params.mean, params.nRows, params.nCols, params.max_rank);
  vector::FVector<Params_t*> pars;
  pars.values = new Params_t*[mus.size()];
  pars.size = mus.size();
  params.step_size = step_sizes[0];
  params.mu = mus[0];
  pars.values[0] = &params;
  for (unsigned i = 1; i < pars.size; i++) {
    Params_t *temp = new Params_t(params);
    temp->step_size = step_sizes[i];
    temp->mu = mus[i];
    pars.values[i] = temp;
  }

  AutoTuner bb (model, pars, nSplits);


  MemoryScan<types::Entry> tscan(test_examps);

  fflush(stdout);
  if (bFileScan) {
    scan::BinaryFileScanner scan(szExampleFile);
    printf("Making scan...\n");
    hogwild::FileScan<scan::BinaryFileScanner, types::Entry>
        fscan(scan, bFileScan);
    fscan.Init();
    printf("Running experiment...\n");
    bb.RunExperiment(nEpochs, wall_clock, fscan, tscan);
  } else {
    MemoryScan<types::Entry> mscan(train_examps);
    bb.RunExperiment(nEpochs, wall_clock, mscan, tscan);
  }

  fflush(stdout);
  if (outputTestFile != NULL) {
    model.OutputToFile(outputTestFile, params);
  }

  fflush(stdout);
  return 0;

}
