

#ifndef HAZY_HOGWILD_BBALL_INL_H
#define HAZY_HOGWILD_BBALL_INL_H



namespace hazy {
namespace hogwild {

template <class Model, class Params, class Example>
class BBall {
 public:
  
  struct BModel {
    Hogwild<Model, Params> *hogs;
    Scan<Example> *scan;
    double *results;
  }

  struct BParams {
  };

  /*! \brief Execute something for each ball
   */
  template <class Exec>
  struct ForEachExec {
    void Execute(BModel &m, const BParams &p, unsigned tid, unsigned tot) {
      hogs[tid].ForEach<Exec>(*m.scan);
    }
  };

  template <class Exec>
  struct AggForEachExec {
    void Execute(BModel &m, const BParams &p, unsigned tid, unsigned tot) {
      m.results[i] = hogs[tid].AggForEach<Exec>(*m.scan);
    }
  };

  template <class Exec, template <class E> class Scan, class Ex>
  void ForEach(Scan<Ex> &scan);

  template <class Exec, template <class E> class Scan, class Ex>
  typename Exec::Aggregate_t AggForEach(Scan<Ex> &scan);




};

} // namespace hogwild
} // namespace hazy
#endif
