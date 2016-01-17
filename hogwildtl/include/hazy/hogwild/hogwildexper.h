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


#ifndef HAZY_HOGWILD_HOGWILDexper_H
#define HAZY_HOGWILD_HOGWILDexper_H

#include "hazy/hogwild/hogwild2-inl.h"

namespace hazy {
namespace hogwild {

/*! A simple class to run timing experiments for ML and analytic tasks
 * \tparam Exper Must have paritcular functions, see src/test/test_hogwild.cc
 */
template <class Exper>
class HogwildExper {
  typedef typename Exper::Model_t Model_t;
  typedef typename Exper::Params_t Params_t;
  typedef typename Exper::Example_t Example_t;

  /*! \brief Wrapper for model training, passed to Hogwild */
  struct TrainExec { 
    static void Execute(Model_t &m, Params_t const &p, Example_t const &e) {
      Exper::TrainModel(m, p, e);
    }
  };

  /*! \brief Wrapper for RMSE computation, passed to Hogwild
   * Implemented the aggregation functions to support RMSE computation
   */
  struct RMSEExec { 
    //! \brief Required (by name) declaration of the type of aggregate
    typedef double Aggregate_t;
    //! \brief Required (by name) declaration of the initial value
    static Aggregate_t Aggregate_Identity;
    static Aggregate_t Execute(Model_t const &m, Params_t const &p, Example_t const &e) {
      return Exper::ComputeLoss(m, p, e);
    }
    static Aggregate_t Merge(Aggregate_t a, Aggregate_t b) { return a + b; }
    static Aggregate_t Aggregate(Aggregate_t a, size_t total) { 
      return std::sqrt(a) / std::sqrt(total); 
    }
  };

 public:
  HogwildExper(Exper &exper, unsigned nthreads) : exper_(exper), tpool_(nthreads), 
      hogwild_(exper.Model(), exper.Params(), tpool_) {
      tpool_.Init();
    }

  void TrainModel();

  double ComputeTestRMSE();
  double ComputeTrainRMSE();

  void RunExperimentWithTest(int epochs);

 private:
  Exper &exper_;
  thread::ThreadPool tpool_;
  Hogwild<Model_t, Params_t> hogwild_;

};

template <class Exper>
  typename HogwildExper<Exper>::RMSEExec::Aggregate_t
  HogwildExper<Exper>::RMSEExec::Aggregate_Identity = 0;

} // namespace hogwild
} // namespace hazy

#endif
