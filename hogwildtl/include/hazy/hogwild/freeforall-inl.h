// Copyright 2012 Chris Re, Victor Bittorf
//
 //Licensed under the Apache License, Version 2.0 (the "License");
 //you may not use this file except in compliance with the License.
 //You may obtain a copy of the License at
 //    http://www.apache.org/licenses/LICENSE-2.0
 //Unless required by applicable law or agreed to in writing, software
 //distributed under the License is distributed on an "AS IS" BASIS,
 //WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 //See the License for the specific language governing permissions and
 //limitations under the License.

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)


#ifndef HAZY_HOGWILD_FREEFORALL_INL_H
#define HAZY_HOGWILD_FREEFORALL_INL_H

#include "hazy/hogwild/tools-inl.h"
#include "hazy/hogwild/hogwild_task.h"

// See for documentation
#include "hazy/hogwild/freeforall.h"

namespace hazy {
namespace hogwild {


template <class Model, class Params, class Exec>
class GoWildDeleg {
 public:
  struct Pack {
    Pack(Model &mm, Params &pp) : m(mm), p(pp) { }
    Model &m;
    Params &p;
  };
  static void GoWildHook(Pack &pack, unsigned tid, unsigned tot) {
    Exec::Execute(pack.m, pack.p, tid, tot);
  }
};


template <class Model, class Params, class Ex, class Exec>
class ForEachDeleg {
 public:
 struct Pack {
    Pack(Model &mm, Params &pp, vector::FVector<Ex> v,
         vector::FVector<size_t> per) 
        : model(mm), params(pp), examps(v), perm(per) { }
    Model &model;
    Params &params;
    vector::FVector<Ex> examps;
    vector::FVector<size_t> perm;
  };

  static void ForEachHook(Pack &pack, unsigned tid, unsigned tot) {
    size_t start = GetStartIndex(pack.examps.size, tid, tot);
    size_t end = GetEndIndex(pack.examps.size, tid, tot);

    size_t *perm = pack.perm.values;
    const Ex *examps = pack.examps.values;
    const Params &params = pack.params;
    Model &model = pack.model;

    for (size_t i = start; i < end; i++) {
      size_t indirect = perm[i];
      Exec::Execute(model, params, examps[indirect]);
    }
  }


  static void ForEach(Pack &task, hazy::thread::ThreadPool &tpool) {
    tpool.Execute(task, ForEachDeleg<Model, Params, Ex, Exec>::ForEachHook);
    tpool.Wait();
  }
};


template <class Model, class Params, class Ex, class Exec>
class AggDelegate {
  typedef typename Exec::Aggregate_t Aggregate_t;
 public:
 struct Pack {
    Pack(Model &mm, Params &pp, vector::FVector<Ex> v,
         vector::FVector<size_t> per) 
        : model(mm), params(pp), examps(v), perm(per) { }
    Model &model;
    Params &params;
    vector::FVector<Ex> examps;
    vector::FVector<size_t> perm;
    
    Aggregate_t *results;
  };

  static void ForEachHook(Pack &pack, unsigned tid, unsigned tot) {
    size_t start = GetStartIndex(pack.examps.size, tid, tot);
    size_t end = GetEndIndex(pack.examps.size, tid, tot);

    size_t *perm = pack.perm.values;
    const Ex *examps = pack.examps.values;
    const Params &params = pack.params;
    Model &model = pack.model;

    Aggregate_t agg = pack.results[tid];

    for (size_t i = start; i < end; i++) {
      size_t indirect = perm[i];
      Aggregate_t agg_temp = Exec::Execute(model, params, examps[indirect]);
      agg = Exec::Merge(agg, agg_temp);
    }
    pack.results[tid] = agg;
  }

  static void BlockHook(Pack &pack, unsigned tid, unsigned tot) {
    const Params &params = pack.params;
    Model &model = pack.model;

    pack.results[tid] = Exec::Execute(model, params, pack.examps, pack.perm,
                                      tid, tot);
  }

  static Aggregate_t ForEach(Pack &pack, hazy::thread::ThreadPool &tpool) {
    // Set up tasks and counters
    Aggregate_t result[tpool.ThreadCount()];
    pack.results = result;

    for (size_t i = 0; i < tpool.ThreadCount(); i++) {
      pack.results[i] = Exec::Aggregate_Identity;
    }

    // Use tasks to populate the aggergate results
    tpool.Execute(pack, AggDelegate<Model, Params, Ex, Exec>::ForEachHook);
    tpool.Wait();

    // Merge the aggregate resgults
    Aggregate_t agg = Exec::Aggregate_Identity;

    for (size_t i = 0; i < tpool.ThreadCount(); i++) {
      agg = Exec::Merge(agg, pack.results[i]);
    }
    return agg;
  }

  static Aggregate_t AggBlockScan(Pack &pack, 
                                  hazy::thread::ThreadPool &tpool) {

    // Set up tasks and counters
    Aggregate_t result[tpool.ThreadCount()];
    pack.results = result;

    for (size_t i = 0; i < tpool.ThreadCount(); i++) {
      pack.results[i] = Exec::Aggregate_Identity;
    }

    // Use tasks to populate the aggergate results
    tpool.Execute(pack, AggDelegate<Model, Params, Ex, Exec>::BlockHook);
    tpool.Wait();

    // Merge the aggregate resgults
    Aggregate_t agg = Exec::Aggregate_Identity;

    for (size_t i = 0; i < tpool.ThreadCount(); i++) {
      agg = Exec::Merge(agg, pack.results[i]);
    }
    return agg;
  }
};


template <class HogwildTask_t>
struct FFATask {
  HogwildTask_t *task;
  vector::FVector<double> *result;
  double (*hook)(HogwildTask_t&, unsigned, unsigned);
};

template <class HogwildTask_t>
void FFADelegate(FFATask<HogwildTask_t> &task, unsigned tid, unsigned tot) {
  assert(task.result->size > tid);
  double res = task.hook(*task.task, tid, tot);
  task.result->values[tid] += res;
}

template <class HogwildTask_t>
void FreeForAll(HogwildTask_t &task, hazy::thread::ThreadPool &tpool,
                double (*hook)(HogwildTask_t&, unsigned, unsigned),
                vector::FVector<double> &result) {
  if (result.size < tpool.ThreadCount()) {
    printf("#threads = %u, size = %lu\n", tpool.ThreadCount(), result.size);
  }
  assert(result.size >= tpool.ThreadCount());

  FFATask<HogwildTask_t> tsk;
  tsk.task = &task;
  tsk.result = &result;
  tsk.hook = hook;

  tpool.Execute(tsk, FFADelegate<HogwildTask_t>);
  tpool.Wait();
}

template <class Model, class Params, class Example, class Scan>
size_t FFAScan(Model &m, Params &p, Scan &scan, 
               hazy::thread::ThreadPool &tpool,
     double (*hook)(HogwildTask<Model, Params, Example>&, unsigned, unsigned),
               vector::FVector<double> &result) {

  size_t count = 0;

  HogwildTask<Model, Params, Example> task;
  task.model = &m;
  task.params = &p;

  while (scan.HasNext()) {
    ExampleBlock<Example> &ex = scan.Next();
    task.block = &ex;
    count += ex.ex.size;

    FreeForAll(task, tpool, hook, result);
  }
  return count;
}


/*
template <class HogwildTask_t, class ATypeSpace>
void FFAForEach() {
  if (result.size < tpool.ThreadCount()) {
    printf("#threads = %u, size = %lu\n", tpool.ThreadCount(), result.size);
  }
  assert(result.size >= tpool.ThreadCount());

  ForEachDelegate<ATypeSpace>::ForEachDelg::Pack 
      pack(*task.model, *task.params, task.block->ex, task.block->perm);

  tpool.Wait();
}


template <class Model, class Params, class Example, class Scan>
size_t FFAScan(Model &m, Params &p, Scan &scan, 
               hazy::thread::ThreadPool &tpool,
     double (*hook)(HogwildTask<Model, Params, Example>&, unsigned, unsigned),
               vector::FVector<double> &result) {

  size_t count = 0;

  HogwildTask<Model, Params, Example> task;
  task.model = &m;
  task.params = &p;

  while (scan.HasNext()) {
    ExampleBlock<Example> &ex = scan.Next();
    task.block = &ex;
    count += ex.ex.size;

    FreeForAll(task, tpool, hook, result);
  }
  return count;
}

};


template <class HogwildTask_t>
void FFADelegate(FFATask<HogwildTask_t> &task, unsigned tid, unsigned tot) {
  assert(task.result->size > tid);
  double res = task.hook(*task.task, tid, tot);
  task.result->values[tid] += res;
}
template <class HogwildTask_t>
void FreeForAll(HogwildTask_t &task, hazy::thread::ThreadPool &tpool,
                double (*hook)(HogwildTask_t&, unsigned, unsigned),
                vector::FVector<double> &result) {
  if (result.size < tpool.ThreadCount()) {
    printf("#threads = %u, size = %lu\n", tpool.ThreadCount(), result.size);
  }
  assert(result.size >= tpool.ThreadCount());

  FFATask<HogwildTask_t> tsk;
  tsk.task = &task;
  tsk.result = &result;
  tsk.hook = hook;

  tpool.Execute(tsk, FFADelegate<HogwildTask_t>);
  tpool.Wait();
}

template <class Model, class Params, class Example, class Scan>
size_t FFAScan(Model &m, Params &p, Scan &scan, 
               hazy::thread::ThreadPool &tpool,
     double (*hook)(HogwildTask<Model, Params, Example>&, unsigned, unsigned),
               vector::FVector<double> &result) {

  size_t count = 0;

  HogwildTask<Model, Params, Example> task;
  task.model = &m;
  task.params = &p;

  while (scan.HasNext()) {
    ExampleBlock<Example> &ex = scan.Next();
    task.block = &ex;
    count += ex.ex.size;

    FreeForAll(task, tpool, hook, result);
  }
  return count;
}
*/

} // namespace hogwild
} // namespace hazy
#endif
