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

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_THREAD_THREAD_POOL_INL_H
#define HAZY_THREAD_THREAD_POOL_INL_H

#include <assert.h>
#include <sched.h>
#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <streambuf>
#include <algorithm>

// See for documentation
#include "thread_pool.h"

namespace hazy {
namespace thread {

namespace __threadpool {
// Hook which pthread_create will call for us
// This simply calls back into the thread pool
void* RunThread(void* thread_meta) {
  ThreadMeta &meta = *static_cast<ThreadMeta*>(thread_meta);
  ThreadPool *tp = static_cast<ThreadPool*>(meta.tpool);
  tp->ThreadLoop(meta);
  return NULL;
}

// A way to keep the types for the caller, this will be the call back
// in the thread pool (callback_)
template <class T>
void Invoke(ThreadPool &tp, unsigned thread_id) {
  tp.ThreadCallback<T>(thread_id);
}
} // namespace __threadpool

ThreadPool::~ThreadPool() {
  if (threads_ == NULL) {
    return;
  }
  delete [] threads_;
  delete [] metas_; // FIXME: this might be deleted before threads exit
  delete [] thread_core_mapping_;
  delete [] thread_node_mapping_;
  delete [] thread_phycore_mapping_;
  delete [] cpuids_;
  barrier_destroy(&ready_);
  barrier_destroy(&finished_);
}

void ThreadPool::Init() {
  hook_ = NULL;
  callback_ = NULL;
  arg_ = NULL;

  // +1 because 'main' thread will join the barriers
  barrier_init(&ready_, NULL, n_threads_+1);
  barrier_init(&finished_, NULL, n_threads_+1);

  metas_ = new ThreadMeta[n_threads_];

  SetExitFlags(false);
  for (unsigned i = 0; i < n_threads_; i++) {
    metas_[i].thread_id = i;
    metas_[i].ready = &ready_;
    metas_[i].finished = &finished_;
    metas_[i].tpool = this;
    metas_[i].binded = false;
  }

  if(numa_available() < 0) {
    printf("System does not support NUMA API!\n");
    exit(0);
  }
  ncpus_ = numa_num_task_cpus();
  nnodes_ = numa_max_node() + 1;
  nphycpus_ = 0;
  printf("We are running on %d nodes and %d CPUs\n", nnodes_, ncpus_);
  cpuids_ = new std::vector<std::vector<int> >[nnodes_];
  GetTopology();
  printf("CPU Topology: \n");
  for (unsigned i = 0; i < nnodes_; ++i) {
    printf("node %d:\t", i);
    for (std::vector<std::vector<int> >::const_iterator cpu = cpuids_[i].begin(); cpu != cpuids_[i].end(); ++cpu) {
      printf("[ ");
      for (std::vector<int>::const_iterator t = (*cpu).begin(); t != (*cpu).end(); ++t)
        printf("%d ", *t);
      printf("] ");
    }
    putchar('\n');
    nphycpus_ += cpuids_[i].size();
  }
  printf("%d physical cores total.\n", nphycpus_);
  ConfigThreadAffinity(); 
  threads_ = new pthread_t[n_threads_];
  for (unsigned i = 0; i < n_threads_; i++) {
    pthread_create(&threads_[i], NULL, __threadpool::RunThread,
                   static_cast<void*>(&metas_[i]));
  }
  ready_flag_ = true;
}

void ThreadPool::GetTopology() {
  std::vector<int> known_siblings;
  for (unsigned cpu = 0; cpu < ncpus_; ++cpu) {
    // skip a core if it is a hyper-threaded logical core
    if (std::find(known_siblings.begin(), known_siblings.end(), cpu) != known_siblings.end())
      continue;
    // if this core is not known as a sibling, add it to the core vector of its node
    int node = numa_node_of_cpu(cpu);
    std::vector<int> phy_core;
    // find out the siblings of this core
    std::stringstream path;
    path << "/sys/devices/system/cpu/cpu" << cpu << "/topology/thread_siblings_list";
    std::ifstream f(path.str().c_str());
    if (f) {
      std::string siblings((std::istreambuf_iterator<char>(f)),
		     std::istreambuf_iterator<char>());
      f.close();
      // std::cout << "CPU " << cpu << " siblings:" << siblings << std::endl;
      std::istringstream ss(siblings);
      std::string coreid;
      while(std::getline(ss, coreid, ',')) {
	known_siblings.push_back(atoi(coreid.c_str()));
        phy_core.push_back(atoi(coreid.c_str()));
      }
      /*
      printf("Known siblings when processing core %d\n", cpu); 
      for (std::vector<int>::const_iterator t = known_siblings.begin(); t != known_siblings.end(); ++t)
	printf("%d ", *t);
      putchar('\n');
      */
    }
    else {
      phy_core.push_back(cpu);
    }
    cpuids_[node].push_back(phy_core);
  }
}

void ThreadPool::ConfigThreadAffinity() {
  thread_core_mapping_ = new int[n_threads_];
  thread_node_mapping_ = new int[n_threads_];
  thread_phycore_mapping_ = new int[n_threads_];
  int node_id, core_id, phycore_id;
  last_used_node_ = 0;
  for (unsigned i = 0; i < n_threads_; ++i) {
    AssignThreadAffinity(i, &node_id, &core_id, &phycore_id);
    thread_core_mapping_[i] = core_id;
    thread_node_mapping_[i] = node_id;
    last_used_node_ = std::max(last_used_node_, (unsigned)node_id);
    thread_phycore_mapping_[i] = phycore_id;
    printf("Thread %d mapped to core %d (phycore %d) on node %d\n", i, core_id, phycore_id, node_id);
  }  
  last_used_node_ += 1;
}

int ThreadPool::GetThreadCoreAffinity(unsigned thread_id) const {
  if (thread_core_mapping_ != NULL && thread_id < n_threads_)
    return thread_core_mapping_[thread_id];
  else
    return -1;
}

int ThreadPool::GetThreadPhyCoreAffinity(unsigned thread_id) const {
  if (thread_phycore_mapping_ != NULL && thread_id < n_threads_)
    return thread_phycore_mapping_[thread_id];
  else
    return -1;
}

int ThreadPool::GetThreadNodeAffinity(unsigned thread_id) const {
  if (thread_node_mapping_ != NULL && thread_id < n_threads_)
    return thread_node_mapping_[thread_id];
  else
    return -1;
}

void ThreadPool::BindToCPU(ThreadMeta &meta) {
  struct bitmask * cpu_mask;
  int cpu = thread_core_mapping_[meta.thread_id];
  // printf("Binding thread %d to CPU %d\n", meta.thread_id, cpu);
  cpu_mask = numa_allocate_cpumask();
  numa_bitmask_setbit(cpu_mask, cpu);
  numa_sched_setaffinity(0, cpu_mask);
  numa_free_cpumask(cpu_mask);
}

void ThreadPool::AssignThreadAffinity(unsigned thread_id, int * node_id, int * core_id, int * phycore_id) {
  int node = -1;
  int ht_id = thread_id / nphycpus_;
  int coreid = thread_id % nphycpus_; 
  int i;
  for (i = 0; i <= coreid; i += cpuids_[++node].size());
  i = coreid - i + cpuids_[node].size();
  // printf("i = %d, node = %d, ht_id = %d, coreid = %d\n", i, node, ht_id, coreid);
  int core = cpuids_[node][i][ht_id % cpuids_[node][i].size()];
  int phycore = cpuids_[node][i][0];
  *node_id = node;
  *core_id = core;
  *phycore_id = phycore;
}

void ThreadPool::ThreadLoop(ThreadMeta &meta) {
  while (true) {
    barrier_wait(meta.ready);
    // we want to bind our thread to a specified CPU
    if (!meta.binded) {
      BindToCPU(meta);
      meta.binded = true;
    }
    if (meta.exit_flag) {
      break;
    }
    callback_(*this, meta.thread_id);
    barrier_wait(meta.finished);
  }
}

template <class Task>
void ThreadPool::ThreadCallback(unsigned thread_id) {
  Task *t = static_cast<Task*>(arg_);
  void (*hook)(Task&, unsigned, unsigned) = reinterpret_cast<
      void (*)(Task&, unsigned, unsigned)>(hook_);
  // int cpu;
  // cpu = sched_getcpu();
  // printf("This is thread %u, running on CPU %d\n", thread_id, cpu);
  hook(*t, thread_id, n_threads_);
  // cpu = sched_getcpu();
  // printf("Thread %u exiting hook, running on CPU %d\n", thread_id, cpu);
}

template <class Task>
void ThreadPool::Execute(Task &task, void (*hook)(Task&, unsigned, 
                         unsigned)) {
  assert(ready_flag_);
  // assign each 
  arg_ = &task; 
  hook_ = reinterpret_cast<void*>(hook);
  callback_ = &__threadpool::Invoke<Task>;
  ready_flag_ = false;
  barrier_wait(&ready_);
}

void ThreadPool::Wait() {
  assert(!ready_flag_);
  barrier_wait(&finished_);
  ready_flag_ = true;
}

void ThreadPool::Join() {
  assert(ready_flag_);

  SetExitFlags(true);
  barrier_wait(&ready_);
  for (unsigned i = 0; i < n_threads_; i++) {
    pthread_join(threads_[i], NULL);
  }
}

} // namespace thread
} // namespace hazy
#endif

