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

#ifndef HAZY_HOGWILD_HOGWILD2_INL_H
#define HAZY_HOGWILD_HOGWILD2_INL_H

#include <cmath>
#include <cstdio>

// Contains most of the helper functions for this class
#include "hazy/hogwild/freeforall-inl.h"

// See for documentation
#include "hazy/hogwild/hogwild2.h"

namespace hazy {
namespace hogwild {

template <class Model, class Params>
template <class Exec>
void Hogwild<Model, Params>::WildExecute() {
  typename GoWildDeleg<Model, Params, Exec>::GoWildDelg::Pack 
      pack(model_, params_);
  tpool_.Execute(pack, 
                 GoWildDeleg<Model, Params, Exec>::GoWildDelg::GoWildHook);
  tpool_.Wait();
}

template <class Model, class Params>
template <class Exec, template <class E> class Scan, class Ex>
typename Exec::Aggregate_t Hogwild<Model, Params>::AggForEach( Scan<Ex> &scan) {
  size_t count = 0;
  typename Exec::Aggregate_t result = Exec::Aggregate_Identity;

  // Scan over all of the blocks of examples
  while (scan.HasNext()) {
    ExampleBlock<Ex> &ex = scan.Next();
    count += ex.ex.size;
    
    // Process the block and aggregate the result
    typename AggDelegate<Model, Params, Ex, Exec>::Pack 
        pack(model_, params_, ex.ex, ex.perm);
    // This will invoke Exec::Execute on each individual example
    // Will block until all threads are done
    typename Exec::Aggregate_t part = 
        AggDelegate<Model, Params, Ex, Exec>::ForEach(pack, tpool_);

    result = Exec::Merge(result, part);
  }
  return Exec::Aggregate(result, count);
}

template <class Model, class Params>
template <class Exec, template <class E> class Scan, class Ex>
typename Exec::Aggregate_t Hogwild<Model, Params>::AggBlockScan(
    Scan<Ex> &scan) {
  size_t count = 0;
  typename Exec::Aggregate_t result = Exec::Aggregate_Identity;

  // Scan over all of the blocks of examples
  while (scan.HasNext()) {
    ExampleBlock<Ex> &ex = scan.Next();
    count += ex.ex.size;

    // Process the block and aggregate the result
    typename AggDelegate<Model, Params, Ex, Exec>::Pack 
        pack(model_, params_, ex.ex, ex.perm);
    // This will invoke Exec::Execute on each block of examples
    // Will block until all threads are done
    typename Exec::Aggregate_t part = 
        AggDelegate<Model, Params, Ex, Exec>::AggBlockScan(pack, tpool_);

    result = Exec::Merge(result, part);
  }
  return Exec::Aggregate(result, count);
}

template <class Model, class Params>
template <class Exec, template <class E> class Scan, class Ex>
void Hogwild<Model, Params>::ForEach(Scan<Ex> &scan) {
  // Scan over all of the blocks of examples
  while (scan.HasNext()) {
    ExampleBlock<Ex> &ex = scan.Next();
    typename ForEachDeleg<Model, Params, Ex, Exec>::Pack 
        pack(model_, params_, ex.ex, ex.perm);
    // Invokes Exec::Execute on each individual example
    // Will block until all threads have finished
    ForEachDeleg<Model, Params, Ex, Exec>::ForEach(pack, tpool_);
  }
}

} // namespace hogwild
} // namespace hazy

#endif
