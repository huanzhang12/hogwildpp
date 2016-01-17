#include "gtest/gtest.h"

#include "fake_scan-inl.h"
#include "hazy/hogwild/hogwild2-inl.h"


typedef hazy::hogwild::HogwildTask<int*, int, int> TestTask;

class UpdateTest : public ::testing::TestWithParam<int> {
};

class RMSETest: public ::testing::TestWithParam<int> {
};


class Tester {
 public:
  static double UpdateModel(TestTask &task, unsigned tid, unsigned tot) {
    (*task.model)[tid] = *task.params * tid;
    return 0.0;
  }
  static double TestModel(TestTask &task, unsigned tid, unsigned tot) {
    return *task.params * (tid+1);
    return 0.0;
  }
  static void PostUpate(int &m, int &p) {
  }
};


TEST_P(UpdateTest, CheckThreading) {
  int n = GetParam();

  int model_[n];
  int *model = model_;
  int param = n;

  hazy::hogwild::FakeScan<int> fs(n, 1);
  hazy::thread::ThreadPool tp(n);

  tp.Init();

  hazy::hogwild::Hogwild<int*, int, Tester> hw(model, param, tp);
  hw.UpdateModel(fs);
  for (int i = 0; i < n; i++)
    ASSERT_EQ(i * param, model[i]);

  tp.Join();

}

TEST_P(RMSETest, CheckThreading) {
  int n = GetParam();

  int model_[n];
  int *model = model_;
  int param = n;

  size_t n_pages = 5;
  size_t examps_per_page = 10;
  hazy::hogwild::FakeScan<int> fs(n_pages, examps_per_page);
  hazy::thread::ThreadPool tp(n);

  tp.Init();

  hazy::hogwild::Hogwild<int*, int, Tester> hw(model, param, tp);
  double expected = 0;
  for (int i = 0; i < n; i++)
    expected += (param * (i+1)) * n_pages;
  expected = std::sqrt(expected) / std::sqrt(n_pages * examps_per_page);
  double rmse = hw.ComputeRMSE(fs);
  ASSERT_EQ(expected, rmse);
  tp.Join();
}


INSTANTIATE_TEST_CASE_P(ModelUpdate, UpdateTest, ::testing::Range(1, 4));
INSTANTIATE_TEST_CASE_P(CalculateRMSE, RMSETest, ::testing::Range(1, 4));
