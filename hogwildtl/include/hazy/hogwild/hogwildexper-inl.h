
#ifndef HAZY_HOGWILD_HOGWILDexper_INL_H
#define HAZY_HOGWILD_HOGWILDexper_INL_H

#include "hazy/util/clock.h"

#include "hogwildexper.h"

namespace hazy {
namespace hogwild {

template <class Exper>
void HogwildExper<Exper>::TrainModel() {
  hogwild_.template ForEach<TrainExec>(exper_.TrainScan());
}

template <class Exper>
double HogwildExper<Exper>::ComputeTestRMSE() {
  return hogwild_.template AggForEach<RMSEExec>(exper_.TestScan());
}

template <class Exper>
double HogwildExper<Exper>::ComputeTrainRMSE() {
  return hogwild_.template AggForEach<RMSEExec>(exper_.TrainScan());
}

template <class Exper>
void HogwildExper<Exper>::RunExperimentWithTest(int epochs) {
  util::Clock total;
  util::Clock train;
  util::Clock test;
  util::Clock epoch;

  total.Start();
  double epoch_time = 0;
  for (int e = 0; e < epochs; e++) {
    train.Start(); epoch.Start();
    TrainModel();
    epoch_time = epoch.Stop(); train.Pause();
    test.Start();
    Exper::PostTrain(exper_.Model(), exper_.Params());
    double train_rmse = ComputeTrainRMSE();
    double test_rmse = ComputeTestRMSE();
    test.Pause();
    printf("epoch %03d: total_time: %06.5f epoch_time: %06.5f train_time: %06.5f test_time: %06.5f train_rmse: %e test_rmse: %e\n", e+1, total.Read(), epoch_time, train.Read(), test.Read(), train_rmse, test_rmse);
    Exper::PostTest(exper_.Model(), exper_.Params());
  }
}

} // namespace hogwild
} // namespace hazy

#endif
