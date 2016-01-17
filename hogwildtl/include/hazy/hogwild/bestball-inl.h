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

#ifndef HAZY_HOGWILD_BESTBALL_INL_H
#define HAZY_HOGWILD_BESTBALL_INL_H

#ifdef __APPLE__
#define isnan std::isnan
#endif


#include <cmath>

#include "hazy/vector/fvector.h"
#include "hazy/vector/operations-inl.h"
#include "hazy/thread/thread_pool-inl.h"
#include "hazy/hogwild/hogwild_task.h"
#include "hazy/hogwild/echo_scan-inl.h"
#include "hazy/hogwild/hogwild-inl.h"
#include "hazy/hogwild/freeforall-inl.h"

namespace hazy {
namespace hogwild {

/* \brief Use the `best ball' method to train a model.
 * This is analogous to golf where mutliply 'balls' are hit are everyone
 * plays from the `best ball.' This is a greedy approach which trains
 * several models concurrently (with different parameters) and at the begining
 * of each new epoch starts all to be the best model from the previous epoch.
 * The notion of 'best model' is defined to be the model with the lowest 
 * hamronic mean of the train and test RMSE.
 */
template <class Model, class Params, class Example, class Exec>
class BestBall {
 public:
  /*! \breif The type of Hogwild that will be used for each ball */
  typedef Hogwild<Model, Params, Exec> Hogwild_t;

  /*! \breif The model for each ball in the best ball method
   * This is structured to mirror the normal Hogwild! interface
   */
  struct BBModel {
    vector::FVector<Model*> models;

    BBModel(Model &m, int count) {
      models.size = (unsigned) count;
      models.values = new Model*[count];
      models.values[0] = &m;
      for (int i = 1; i < count; i++) {
        models.values[i] = m.Clone();
      }
    }
  };

  /*! \breif Parameters for each ball in the best ball method
   * This is structured to mirror the normal Hogwild! interface
   */
  struct BBParams {
    vector::FVector<Params*> &params;
    Hogwild_t** hogwilds;
    thread::ThreadPool **tpools;

    BBParams(BBModel &models, vector::FVector<Params*> &ps, unsigned threads) 
          : params(ps) {
      tpools = new thread::ThreadPool*[ps.size];
      hogwilds = new Hogwild_t*[ps.size];
      for (unsigned i = 0; i < ps.size; i++) {
        tpools[i] = new thread::ThreadPool(threads);
        tpools[i]->Init();
        hogwilds[i] = new Hogwild_t(*models.models.values[i], *ps.values[i],
                                    *tpools[i]);
      }
    }
  };

  typedef HogwildTask<BBModel, BBParams, Example> BallTask;

  /*! Worker thread which updates a model based on the thead id
   */
  static double UpdateModel(BallTask &task, unsigned tid, unsigned tot) {
    EchoScan<Example> echo(*task.block);
    task.params->hogwilds[tid]->UpdateModel(echo);
    return 0.0;
  }

  /*! Worker thread which computes the RMSE of a model based on the thead id
   */
  static double TestModel(BallTask &task, unsigned tid, unsigned tot) {
    EchoScan<Example> echo(*task.block);
    double rmse = task.params->hogwilds[tid]->ComputeRMSE(echo);
    double se = rmse * std::sqrt(task.block->ex.size);
    return se * se;
  }

  BestBall(Model &model, vector::FVector<Params*> &ps, unsigned threads) {
    tpool_ = new thread::ThreadPool(ps.size);
    tpool_->Init();
    models_ = new BBModel(model, ps.size);
    params_ = new BBParams(*models_, ps, threads);
    res_.size = ps.size;
    res_.values = new double[res_.size];
    train_.size = ps.size;
    train_.values = new double[res_.size];
    test_.size = ps.size;
    test_.values = new double[res_.size];
  }

  /*! \brief Run a best ball experiment and print the output
   * \param nepochs Number of epcohs to run for
   * \param wall_clock a clock used to report time, started at the start start
   *        of the program.
   * \param trscan the train example scanner
   * \param tescan the test example scanner
   */
  template <class TrainScan, class TestScan>
  void RunExperiment(int nepochs, hazy::util::Clock &wall_clock, 
                     TrainScan &trscan, TestScan &tescan) {
    printf("wall_clock: %.5f    Going Hogwild!\n", wall_clock.Read());
    for (int e = 1; e <= nepochs; e++) {
      UpdateModels(trscan);
      ComputeRMSE(trscan, train_);
      ComputeRMSE(tescan, test_);

      unsigned bb = PickBestBall();

      printf("epoch: %d wall_clock: %.5f best_ball: %u harm_mean: %.5f train_rmse: %.5f test_rmse: %.5f\n", 
             e, wall_clock.Read(), bb, res_.values[bb], train_.values[bb],
             test_.values[bb]);
      fflush(stdout);
    }
  }

  /*! \brief Runs the experiment without test examples, only using train data
   * \param nepochs Number of epcohs to run for
   * \param wall_clock a clock used to report time, started at the start start
   *        of the program.
   * \param trscan the train example scanner
   */
  template <class TrainScan>
  void RunExperiment(int nepochs, hazy::util::Clock &wall_clock, 
                     TrainScan &trscan) {
    printf("wall_clock: %.5f    Going Hogwild!\n", wall_clock.Read());
    for (int e = 1; e <= nepochs; e++) {
      UpdateModels(trscan);
      ComputeRMSE(trscan, train_);
      vector::CopyInto(train_, test_);

      unsigned bb = PickBestBall();

      printf("epoch: %d wall_clock: %.5f best_ball: %u train_rmse: %.5f\n", 
             e, wall_clock.Read(), bb, train_.values[bb]);
      fflush(stdout);
    }
  }


 private:
  BBModel *models_;
  BBParams *params_;
  thread::ThreadPool *tpool_;
  vector::FVector<double> res_; //!< results (for computing RMSE) size=nthreads
  vector::FVector<double> train_; //!< results (for computing RMSE) size=nthreads
  vector::FVector<double> test_; //!< results (for computing RMSE) size=nthreads

  template <class Scan>
  void UpdateModels(Scan &scan) {
    scan.Reset();
    // Discard the result.
    FFAScan(*models_, *params_, scan, *tpool_, UpdateModel, res_);
  }

  /*! \brief Compute the RMSE of each model using the given scanner
   * The results are stored in the given vector, so that the RMSE of
   * model #i is at index i in the vector.
   * \param scan A scanner of examples
   * \param res A vector to place the resulting RMSEs, must already be alloc'd
   */
  template <class Scan>
  void ComputeRMSE(Scan &scan, vector::FVector<double> &res) {
    scan.Reset();
    Zero(res);
    size_t count = FFAScan(*models_, *params_, scan, 
                           *tpool_, TestModel, res);
    for (unsigned i = 0; i < tpool_->ThreadCount(); i++) {
      res.values[i] = std::sqrt(res.values[i]) / std::sqrt(count);
    }
  }
  
  /*! \brief Select the best ball and move all balls to that place
   * This means selecting the ball with the lowest RMSE and copying that
   * model into all of the other models.
   */
  unsigned PickBestBall() {
    // find lowest RMSE
    unsigned lowest = 0;
    for (unsigned i = 0; i < res_.size; i++) {
      res_.values[i] = (2 * train_.values[i] * test_.values[i]) /
          (train_.values[i] + test_.values[i]);
      printf("Ball #%u harmonic mean of RMSEs = %lf\n", i, res_.values[i]);
      if (! isnan(res_.values[i])) {
        lowest = i;
      }
    }
    if (isnan(res_.values[lowest])) {
      printf("All models diverged!!\n");
      assert(false);
    }
    for (unsigned i = 0; i < res_.size; i++) {
      if (isnan(res_.values[lowest])) {
        continue;
      }
      if (res_.values[lowest] > res_.values[i]) {
        lowest = i;
      }
    }

    for (unsigned i = 0; i < res_.size; i++) {
      if (i == lowest) {
        continue;
      }
      models_->models.values[i]->CopyFrom(*models_->models.values[lowest]);
    }
    return lowest;
  }

  void Zero(vector::FVector<double> &res) { for (unsigned i = 0; i < res.size; i++) res.values[i] = 0; }
};

} // namespace hogwild
} // namespace hazy
#endif
