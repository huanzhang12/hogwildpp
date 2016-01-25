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
#include <cstring>
#include <numa.h>

#include "hazy/hogwild/hogwild-inl.h"
#include "hazy/hogwild/numa_memory_scan.h"
#include "hazy/scan/tsvfscan.h"
#include "hazy/scan/binfscan.h"

#include "frontend_util.h"

#include "numasvm/svmmodel.h"
#include "svm/svm_loader.h"
#include "numasvm/svm_exec.h"


// Hazy imports
using namespace hazy;
using namespace hazy::hogwild;
using scan::TSVFileScanner;
using scan::MatlabTSVFileScanner;

using hazy::hogwild::svm::fp_type;


using namespace hazy::hogwild::svm;


template <class Scan>
size_t NumaLoadSVMExamples(Scan &scan, vector::FVector<SVMExample> * nodeex, unsigned nnodes) { 
  size_t nfeats = 0;
#if 0
  for (unsigned i = 0; i < nnodes; ++i) {
    scan.Reset();
    numa_run_on_node(i);
    numa_set_preferred(i);
    nfeats = LoadSVMExamples<Scan>(scan, nodeex[i]);
  }
#else
  numa_run_on_node(0);
  numa_set_preferred(0);
  // The examples on the first node will be loaded from the input file
  nfeats = LoadSVMExamples<Scan>(scan, nodeex[0]);
  // Other nodes need a local copy
  for (unsigned n = 1; n < nnodes; ++n) {
    // scan.Reset();
    numa_run_on_node(n);
    numa_set_preferred(n);
    nodeex[n].size = nodeex[0].size;
    nodeex[n].values = new SVMExample[nodeex[n].size];
    for (size_t i = 0; i < nodeex[n].size; i++) {
      size_t size = nodeex[0].values[i].vector.size;
      nodeex[n].values[i].value = nodeex[0].values[i].value;
      nodeex[n].values[i].vector.size = size;
      nodeex[n].values[i].vector.index = new int[size];
      nodeex[n].values[i].vector.values = new fp_type[size];
      std::memcpy((void *)nodeex[n].values[i].vector.index, nodeex[0].values[i].vector.index, size * sizeof(int));
      std::memcpy((void *)nodeex[n].values[i].vector.values, nodeex[0].values[i].vector.values, size * sizeof(fp_type));
      for (size_t j = 0; j < nodeex[n].values[i].vector.size; j++) {
	// assert(nodeex[n].values[i].vector.index[j] >= 0);
	// assert(nodeex[n].values[i].vector.index[j] <= (long int)nfeats - 1);
      }
    }
  }
#endif
  numa_run_on_node(-1);
  // numa_set_preferred(-1);
  numa_set_localalloc();
  return nfeats;
}

fp_type SolveBeta(int n) {
  fp_type start = 0.6;
  fp_type end = 1.0;
  fp_type mid = 0.5;
  fp_type err = 0;
  if (n >= 2) {
    do {
      mid = (start + end) / 2;
      err = pow(mid, n) + mid - 1;
      if (err > 0) {
	end = mid; 
      }
      else {
	start = mid;
      }
    } while(fabs(err) > 0.001);
  }
  printf("Beta for n=%d is %f (err %f)\n", n, mid, err);
  return mid;
}

void CreateNumaSVMModel(NumaSVMModel * &node_m, size_t nfeats, hazy::thread::ThreadPool &tpool, unsigned nthreads) {
  // Build the weight update chain
  int * thread_to_weights_mapping = new int[nthreads];
  int * next_weights = new int[nthreads];
  /* weight update policy: 
     Hyper-threaded core: update the same set of weights in the same physical core;
     Physical core: has a unique set of weights
  */
  for (unsigned i = 0; i < nthreads; ++i) {
    thread_to_weights_mapping[i] = tpool.GetThreadPhyCoreAffinity(i);
  }
  /* next-weight policy:
     Hyper-threaded core: Do not update next set of weights.
     Physical core: update the weights of next core, only cross node boundary once
  */
  unsigned phycpu_count = tpool.PhyCPUCount();
  unsigned weights_count = nthreads > phycpu_count ? phycpu_count : nthreads;
  for (unsigned i = 0; i < nthreads; ++i) {
    int phy_core = tpool.GetThreadPhyCoreAffinity(i);
    int logical_core = tpool.GetThreadCoreAffinity(i);
    if (phy_core == logical_core) {
      next_weights[i] = (phy_core + 1) % weights_count;
    }
    else {
      next_weights[i] = -1;
    }
  }
  if (nthreads == 1) {
    next_weights[0] = -1;
  }
  /* Now create the Model array 
  */
  numa_run_on_node(0);
  numa_set_preferred(0);
  int * atomic_ptr = new int;
  int atomic_mask = (1 << (sizeof(int) * 8 - (weights_count - 1 ? __builtin_clz(weights_count - 1) : 32))) - 1;
  node_m = new NumaSVMModel[weights_count];
  printf("Model array allocated at %p\n", node_m);
  for (unsigned i = 0; i < weights_count; ++i) {
    numa_run_on_node(numa_node_of_cpu(i));
    numa_set_preferred(numa_node_of_cpu(i));
    printf("Allocating memory for core %d at node %d\n", i, numa_node_of_cpu(i));
    node_m[i].AllocateModel(nfeats);
    node_m[i].atomic_ptr = atomic_ptr;
    node_m[i].atomic_mask = atomic_mask;
    node_m[i].thread_to_weights_mapping = thread_to_weights_mapping;
    node_m[i].next_weights = next_weights;
    if (i == weights_count - 1) {
      node_m[i].atomic_inc_value = atomic_mask - weights_count + 2;
    }
    else {
      node_m[i].atomic_inc_value = 1;
    }
  }
  numa_run_on_node(-1);
  // numa_set_preferred(-1);
  numa_set_localalloc();
}

int main(int argc, char** argv) {
  hazy::util::Clock wall_clock;
  wall_clock.Start();
  //Benchmark::StartExperiment(argc, argv);

  bool matlab_tsv = false;
  bool loadBinary = false;
  unsigned nepochs = 20;
  unsigned nthreads = 1;
  float mu = 1.0, step_size = 5e-2, step_decay = 0.8;
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
        mu = atof(optarg);
        break;
      case 'e':
        nepochs = atoi(optarg);
        break;
      case 'i':
        step_size = atof(optarg);
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

  // we initialize thread pool here because we need CPU topology information
  hazy::thread::ThreadPool tpool(nthreads);
  tpool.Init();
  unsigned nnodes = tpool.NodeCount();
  unsigned phycpu_count = tpool.PhyCPUCount();
  SVMParams tp (step_size, step_decay, mu, SolveBeta(nthreads > phycpu_count ? phycpu_count : nthreads), &tpool);
  
  vector::FVector<SVMExample> * node_train_examps = new vector::FVector<SVMExample>[nnodes];
  vector::FVector<SVMExample> * node_test_examps = new vector::FVector<SVMExample>[nnodes];

  size_t nfeats;

  if (loadBinary) {
    printf("Loading binary file...\n");
    scan::BinaryFileScanner scan(szExampleFile);
    nfeats = NumaLoadSVMExamples(scan, node_train_examps, nnodes);
    printf("Loaded binary file!\n");
  } else if (matlab_tsv) {
    MatlabTSVFileScanner scan(szExampleFile);
    nfeats = NumaLoadSVMExamples(scan, node_train_examps, nnodes);
  } else {
    TSVFileScanner scan(szExampleFile);
    nfeats = NumaLoadSVMExamples(scan, node_train_examps, nnodes);
  }
  if (matlab_tsv) {
    MatlabTSVFileScanner scantest(szTestFile);
    NumaLoadSVMExamples(scantest, node_test_examps, nnodes);
  } else {
    TSVFileScanner scantest(szTestFile);
    NumaLoadSVMExamples(scantest, node_test_examps, nnodes);
  }

  unsigned *degs = new unsigned[nfeats];
  printf("Loaded %lu examples\n", nfeats);
  for (size_t i = 0; i < nfeats; i++) {
    degs[i] = 0;
  }
  CountDegrees(node_train_examps[0], degs);
  tp.degrees = degs;
  tp.ndim = nfeats;

//  hogwild::freeforall::FeedTrainTest(memfeed.GetTrough(), nepochs, nthreads);
  NumaSVMModel * node_m;
  CreateNumaSVMModel(node_m, nfeats, tpool, nthreads); 
  NumaSVMModel m;
  m.AllocateModel(nfeats);
  NumaMemoryScan<SVMExample> mscan(node_train_examps, nnodes);
  Hogwild<NumaSVMModel, SVMParams, NumaSVMExec>  hw(node_m[0], tp, tpool);
  NumaMemoryScan<SVMExample> tscan(node_test_examps, nnodes);

  hw.RunExperiment(nepochs, wall_clock, mscan, tscan);

  return 0;
}

