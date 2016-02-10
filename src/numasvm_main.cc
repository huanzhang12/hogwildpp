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
#include <set>
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
  if (!n)
    mid = .0; 
  return mid;
}


void PrintNumaMemStats() {
  int node_count = numa_max_node() + 1;
  for (int i = 0; i < node_count; ++i) {
    long fr;
    unsigned long sz = numa_node_size(i, &fr);
    printf("Node %d: %13lu Available, %13ld Free\n", i, sz, fr);
  }
}


int CreateNumaPerCoreCentralSVMModel(NumaSVMModel * &node_m, size_t nfeats, hazy::thread::ThreadPool &tpool, unsigned nthreads, int update_delay) {
  // Build the weight update chain
  int * thread_to_weights_mapping = new int[nthreads];
  int * next_weights = new int[nthreads];
  /* weight update policy: 
     Only node 0 has a unique w
  */
  for (unsigned i = 0; i < nthreads; ++i) {
    thread_to_weights_mapping[i] = tpool.GetThreadPhyCoreAffinity(i);
  }
  /* next-weight policy:
     All threads update to the single w
  */
  int node_count = tpool.NodeCount();
  unsigned phycpu_count = tpool.PhyCPUCount();
  int weights_count = nthreads > phycpu_count ? phycpu_count : nthreads;
  for (unsigned i = 0; i < nthreads; ++i) {
    // all threads (except for hyper-threading cores) update the common w
    int phy_core = tpool.GetThreadPhyCoreAffinity(i);
    int logical_core = tpool.GetThreadCoreAffinity(i);
    if (phy_core == logical_core) {
      next_weights[i] = weights_count;
    }
    else {
      next_weights[i] = -1;
    }
  }
  /* Now create the Model array, per-core, central  
  */
  numa_run_on_node(0);
  numa_set_preferred(0);
  int * atomic_ptr = new int ();
  int atomic_mask = (1 << (sizeof(int) * 8 - (weights_count - 1 ? __builtin_clz(weights_count - 1) : 32))) - 1;
  node_m = new NumaSVMModel[weights_count + 1];
  printf("Model array allocated at %p\n", node_m);
  printf("Allocating memory for main vector\n"); 
  node_m[weights_count].AllocateModel(nfeats);
  node_m[weights_count].atomic_ptr = nullptr;
  node_m[weights_count].atomic_mask = 0;
  node_m[weights_count].thread_to_weights_mapping = nullptr;
  node_m[weights_count].next_weights = nullptr;
  PrintNumaMemStats();
  for (int i = 0; i < weights_count; ++i) {
    numa_run_on_node(numa_node_of_cpu(i));
    numa_set_preferred(numa_node_of_cpu(i));
    printf("Allocating memory for core %d on node %d\n",i, numa_node_of_cpu(i)); 
    node_m[i].AllocateModel(nfeats);
    PrintNumaMemStats();
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
    node_m[i].update_atomic_counter = update_delay;
  }
  numa_run_on_node(-1);
  numa_set_localalloc();
  return weights_count;
}

int CreateNumaPerNodeCentralSVMModel(NumaSVMModel * &node_m, size_t nfeats, hazy::thread::ThreadPool &tpool, unsigned nthreads, int update_delay) {
  // Build the weight update chain
  int * thread_to_weights_mapping = new int[nthreads];
  int * next_weights = new int[nthreads];
  /* weight update policy: 
     Each node has a separated node
  */
  for (unsigned i = 0; i < nthreads; ++i) {
    thread_to_weights_mapping[i] = tpool.GetThreadNodeAffinity(i);
  }
  /* next-weight policy:
     The first thread of a node is responsible for communicating with adjacent node
  */
  int node_count = tpool.NodeCount();
  std::set<int> node_set;
  for (unsigned i = 0; i < nthreads; ++i) {
    int node = tpool.GetThreadNodeAffinity(i);
    if (node_set.find(node) == node_set.end()) {
      node_set.insert(node);
      next_weights[i] = -2; // just a flag, will be replace by max. nodes
    }
    else {
      next_weights[i] = -1;
    }
  }
  if (nthreads == 1) {
    next_weights[0] = -1;
  }
  int max_node = *node_set.rbegin() + 1;
  int weights_count = max_node > node_count ? node_count : max_node;
  for (unsigned i = 0; i < nthreads; ++i) {
    if (max_node == 1) 
      next_weights[i] = -1;
   else if (next_weights[i] == -2) {
      // point to the last weights shared by all nodes
      next_weights[i] = weights_count;
   }
  }
  /* Now create the Model array: per-node, central 
  */
  numa_run_on_node(0);
  numa_set_preferred(0);
  int * atomic_ptr = new int ();
  int atomic_mask = (1 << (sizeof(int) * 8 - (weights_count - 1 ? __builtin_clz(weights_count - 1) : 32))) - 1;
  node_m = new NumaSVMModel[weights_count + 1];
  printf("Model array allocated at %p\n", node_m);
  node_m[weights_count].AllocateModel(nfeats);
  node_m[weights_count].atomic_ptr = nullptr;
  node_m[weights_count].atomic_mask = 0;
  node_m[weights_count].thread_to_weights_mapping = nullptr;
  node_m[weights_count].next_weights = nullptr;
  PrintNumaMemStats();
  for (int i = 0; i < weights_count; ++i) {
    numa_run_on_node(i);
    numa_set_preferred(i);
    printf("Allocating memory for node %d\n",i);
    node_m[i].AllocateModel(nfeats);
    PrintNumaMemStats();
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
    node_m[i].update_atomic_counter = update_delay;
  }
  numa_run_on_node(-1);
  numa_set_localalloc();
  return weights_count;
}


int CreateNumaClusterRoundRobinRingSVMModel(NumaSVMModel * &node_m, size_t nfeats, hazy::thread::ThreadPool &tpool, unsigned nthreads, unsigned cluster_size, int update_delay) {
  // Build the weight update chain
  int * thread_to_weights_mapping = new int[nthreads];
  int * next_weights = new int[nthreads];
  unsigned phycpu_count = tpool.PhyCPUCount();
  int weights_count = nthreads > phycpu_count ? phycpu_count : nthreads;
  int node_count = tpool.NodeCount();
  int cluster_count = weights_count / cluster_size;
  assert((nthreads % cluster_size == 0) && "Total number of threads must be a multiple of cluster size\n");
  assert((nthreads > phycpu_count ? (phycpu_count % cluster_size == 0) : 1) && 
         "When total number threads is greater than core count, core count must by a multiple of cluster size\n");
  /* weight update policy: 
     Each physical core has a separated model data structure, but some of them might share the same w
  */
  for (unsigned i = 0; i < nthreads; ++i) {
    thread_to_weights_mapping[i] = ((i % phycpu_count) % cluster_size) * cluster_count + ((i % phycpu_count) / cluster_size);
  }
  /* next-weight policy:
     Threads of a node take turns communicating with adjacent node
     0, 4, 1, 5, 2, 6, 3, 7
     We assume that first N threads will be allocated to physical cores, 
     (N = total #CPU)
  */
  for (unsigned i = 0; i < nthreads; ++i) {
    if (i < phycpu_count) {
      next_weights[i] = (thread_to_weights_mapping[i] + 1) % weights_count;
    }
    else {
      next_weights[i] = -1;
    }
  }

 /* Now create the Model array: per-cluster, ring 
  */
  numa_run_on_node(0);
  numa_set_preferred(0);
  int * atomic_ptr = new int ();
  int atomic_mask = (1 << (sizeof(int) * 8 - (weights_count - 1 ? __builtin_clz(weights_count - 1) : 32))) - 1;
  node_m = new NumaSVMModel[weights_count]; // some of them are just pointers to other weights
  printf("Model array allocated at %p\n", node_m);
  PrintNumaMemStats();
  for (int i = 0; i < weights_count; ++i) {
    int thread_id = (i % cluster_count) * cluster_size + (i / cluster_count);
    int node = tpool.GetThreadNodeAffinity(thread_id);
    numa_run_on_node(node);
    numa_set_preferred(node);
    if (i / cluster_count == 0) {
      printf("Allocating memory for weight %d (thread %d) on node %d\n", i, thread_id, node);
      node_m[i].AllocateModel(nfeats);
      PrintNumaMemStats();
    }
    else {
      node_m[i].MirrorModel(node_m[i % cluster_count]);
    }
    node_m[i].atomic_ptr = atomic_ptr;
    node_m[i].atomic_mask = atomic_mask;
    if (i == weights_count - 1) {
      node_m[i].atomic_inc_value = atomic_mask - weights_count + 2;
    }
    else {
      node_m[i].atomic_inc_value = 1;
    }
    node_m[i].thread_to_weights_mapping = thread_to_weights_mapping;
    node_m[i].next_weights = next_weights;
    node_m[i].update_atomic_counter = update_delay;
  }
  numa_run_on_node(-1);
  numa_set_localalloc();
  return weights_count;
}

int CreateNumaPerNodeRingSVMModel(NumaSVMModel * &node_m, size_t nfeats, hazy::thread::ThreadPool &tpool, unsigned nthreads, int update_delay) {
  // Build the weight update chain
  int * thread_to_weights_mapping = new int[nthreads];
  int * next_weights = new int[nthreads];
  /* weight update policy: 
     Each node has a separated node
  */
  for (unsigned i = 0; i < nthreads; ++i) {
    thread_to_weights_mapping[i] = tpool.GetThreadNodeAffinity(i);
  }
  /* next-weight policy:
     The first thread of a node is responsible for communicating with adjacent node
  */
  int node_count = tpool.NodeCount();
  std::set<int> node_set;
  for (unsigned i = 0; i < nthreads; ++i) {
    int node = tpool.GetThreadNodeAffinity(i);
    if (node_set.find(node) == node_set.end()) {
      node_set.insert(node);
      next_weights[i] = node + 1;
    }
    else {
      next_weights[i] = -1;
    }
  }
  if (nthreads == 1) {
    next_weights[0] = -1;
  }
  int max_node = *node_set.rbegin() + 1;
  int weights_count = max_node > node_count ? node_count : max_node;
  for (unsigned i = 0; i < nthreads; ++i) {
    if (max_node == 1) 
      next_weights[i] = -1;
   else if (next_weights[i] == max_node)
      next_weights[i] = 0;
  }
  /* Now create the Model array: per-node, ring 
  */
  numa_run_on_node(0);
  numa_set_preferred(0);
  int * atomic_ptr = new int ();
  int atomic_mask = (1 << (sizeof(int) * 8 - (weights_count - 1 ? __builtin_clz(weights_count - 1) : 32))) - 1;
  node_m = new NumaSVMModel[weights_count];
  printf("Model array allocated at %p\n", node_m);
  PrintNumaMemStats();
  for (int i = 0; i < weights_count; ++i) {
    numa_run_on_node(i);
    numa_set_preferred(i);
    printf("Allocating memory for node %d\n",i);
    node_m[i].AllocateModel(nfeats);
    PrintNumaMemStats();
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
    node_m[i].update_atomic_counter = update_delay;
  }
  numa_run_on_node(-1);
  numa_set_localalloc();
  return weights_count;
}

int CreateNumaPerCoreRingSVMModel(NumaSVMModel * &node_m, size_t nfeats, hazy::thread::ThreadPool &tpool, unsigned nthreads, int update_delay) {
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
  /* Now create the Model array, per-core, ring
  */
  numa_run_on_node(0);
  numa_set_preferred(0);
  int * atomic_ptr = new int ();
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
    node_m[i].update_atomic_counter = update_delay;
  }
  numa_run_on_node(-1);
  numa_set_localalloc();
  return weights_count;
}

void PrintWeights(NumaSVMModel * node_m, int weights_count, int nthreads, hazy::thread::ThreadPool &tpool) {
  printf("Thread to weights map:\n");
  NumaSVMModel &model = node_m[0];
  for (int i = 0; i < nthreads; ++i) {
    int weights_index = model.thread_to_weights_mapping[i];
    int next_weights = model.next_weights[i];
    NumaSVMModel * const m = &node_m[weights_index];
    NumaSVMModel * const next_m = next_weights >= 0 ? &node_m[next_weights] : NULL;
    int atomic_inc_value = m->atomic_inc_value;
    int atomic_mask = m->atomic_mask;
    printf("Thread %2d (node %2d phycore %2d core %2d): "
	   "%2d->%2d at %p->%p, (atomic+%2d) & %2x\n", 
	   i, tpool.GetThreadNodeAffinity(i), tpool.GetThreadPhyCoreAffinity(i), 
           tpool.GetThreadCoreAffinity(i), weights_index, next_weights,
           m->weights.values, next_weights >= 0 ? next_m->weights.values: NULL,
           atomic_inc_value, atomic_mask);
  }
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
  bool perCore = true;
  bool useRing = true;
  int cluster_size = 0;
  int update_delay = 256;
  double tolerance = 1e-2;
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
    {"percore", required_argument, NULL, 'p', "Using per-core weights instead of per-node (default: on)"},
    {"ring", required_argument, NULL, 'g', "Using the ring update scheme, otherwise use a common w (default: on)"},
    {"update_delay", required_argument, NULL, 't', "Number of iterations before pass the token to the next thread (default: 256)"},
    {"cluster_size", required_argument, NULL, 'c', "Threads in a cluster share the same weights"},
    {"tolerance", required_argument, NULL, 'o', "error tolerance when doing gradient update (default 1e-2)"},
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
      case 'p':
        perCore = atoi(optarg);
        break;
      case 'g':
        useRing = atoi(optarg);
        break;
      case 't':
        update_delay = atoi(optarg);
        break;
      case 'c':
        cluster_size = atoi(optarg);
        break;
      case 'o':
        tolerance = atof(optarg);
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
  
  unsigned nnodes = tpool.UsedNodeCount();
  printf("%d threads will be running on %d nnodes\n", nthreads, nnodes);
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
  if (loadBinary) {
    printf("Loading binary file...\n");
    scan::BinaryFileScanner scantest(szTestFile);
    NumaLoadSVMExamples(scantest, node_test_examps, nnodes);
    printf("Loaded binary file!\n");
  } else if (matlab_tsv) {
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

  NumaSVMModel * node_m;
  int weights_count;
  fp_type beta, lambda;
  if (cluster_size) {
     weights_count = CreateNumaClusterRoundRobinRingSVMModel(node_m, nfeats, tpool, nthreads, cluster_size, update_delay);
     beta = SolveBeta(weights_count / cluster_size);
     lambda = 1 - pow(beta, weights_count / cluster_size - 1);
  }
  else {
    if (useRing) {
      if (perCore) {
	weights_count = CreateNumaPerCoreRingSVMModel(node_m, nfeats, tpool, nthreads, update_delay); 
      }
      else {
	weights_count = CreateNumaPerNodeRingSVMModel(node_m, nfeats, tpool, nthreads, update_delay); 
      }
      beta = SolveBeta(weights_count);
      lambda = 1 - pow(beta, weights_count - 1);
    }
    else {
      if (perCore) {
	weights_count = CreateNumaPerCoreCentralSVMModel(node_m, nfeats, tpool, nthreads, update_delay); 
      }
      else {
	weights_count = CreateNumaPerNodeCentralSVMModel(node_m, nfeats, tpool, nthreads, update_delay); 
      }
      beta = 1.0;
      lambda = 1.0;
    }
  }
  printf("weights_count=%d, beta=%f, lambda=%f\n", weights_count, beta, lambda);
  PrintWeights(node_m, weights_count, nthreads, tpool);
  SVMParams tp (step_size, step_decay, mu, beta, lambda, weights_count, useRing, update_delay, tolerance, &tpool);
  CountDegrees(node_train_examps[0], degs);
  tp.degrees = degs;
  tp.ndim = nfeats;

//  hogwild::freeforall::FeedTrainTest(memfeed.GetTrough(), nepochs, nthreads);
  NumaMemoryScan<SVMExample> mscan(node_train_examps, nnodes);
  Hogwild<NumaSVMModel, SVMParams, NumaSVMExec>  hw(node_m[0], tp, tpool);
  NumaMemoryScan<SVMExample> tscan(node_test_examps, nnodes);

  hw.RunExperiment(nepochs, wall_clock, mscan, tscan);
  
  return 0;
}

