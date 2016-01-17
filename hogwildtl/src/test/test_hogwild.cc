

#include "hazy/hogwild/hogwildexper-inl.h"
#include "hazy/hogwild/memory_scan.h"

using namespace hazy;
using namespace hogwild;
using namespace vector;

class Exper {
 public:
  typedef double Model_t;
  typedef double Params_t;
  typedef double Example_t;
  typedef double Aggregate_t;

  static void TrainModel(double &m, double p, double e) {
    m += p * e;
  }

  static double ComputeLoss(double m, double p, double e) {
    return m - p * e;
  }

  static void PostTrain(double &m, double &p) { 
    printf("Model = %f\n", m); 
  }
  static void PostTest(double &m, double &p) { }

  Exper(double p, FVector<double> v) : v_(v), p_(p), m_(0), scan_(v) { }

  Model_t& Model() { return m_; }
  Params_t& Params() { return p_; }

  MemoryScan<double>& TrainScan() { scan_.Reset(); return scan_; }
  MemoryScan<double>& TestScan() { scan_.Reset(); return scan_; }


  FVector<double> v_;
  double p_;
  double m_;

  MemoryScan<double> scan_;

};

int main() {

  double arr[] = {1, 2, 3, 4};
  FVector<double> fv (arr, 4);

  Exper ex (5, fv);

  HogwildExper<Exper> hwexper(ex, 2);

  hwexper.RunExperimentWithTest(2);
}
