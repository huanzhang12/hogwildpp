
#include <vector>

#include "hazy/vector/fvector.h"
#include "svmmodel.h"

namespace hazy {
namespace hogwild {
namespace svm {

template <class Scan>
size_t LoadSVMExamples(Scan &scan, vector::FVector<SVMExample> &ex) {
  std::vector<SVMExample> examps;
  int lastrow = -1;
  double rating = 0.0;
  std::vector<fp_type> data;
  std::vector<int> index;

  int max_col = 0;

  while (scan.HasNext()) {
    const types::Entry &e = scan.Next();
    if (lastrow == -1) {
      lastrow = e.row;
    }
    if ((lastrow != e.row) || (!scan.HasNext())) {
      // finish off the previous vector and start a new one
      lastrow = e.row;
      fp_type *d = new fp_type[data.size()];
      int *i = new int[data.size()];
      for (size_t j = 0; j < data.size(); j++) {
        d[j] = data[j];
        i[j] = index[j];
      }
      SVMExample temp(rating, d, i, data.size());
      examps.push_back(temp);
      rating = 0.0;
      data.clear();
      index.clear();
    }

    if (e.col < 0) {
      rating = e.rating;
    } else {
      if (e.col > max_col) {
        max_col = e.col;
      }
      data.push_back(e.rating);
      index.push_back(e.col);
    }
  }

  // Copy from temp vector into persistent memory
  ex.size = examps.size();
  ex.values = new SVMExample[ex.size];
  for (size_t i = 0; i < ex.size; i++) {
    new (&ex.values[i]) SVMExample(examps[i]);
    for (size_t j = 0; j < ex.values[i].vector.size; j++) {
      assert(ex.values[i].vector.index[j] >= 0);
      assert(ex.values[i].vector.index[j] <= max_col);
    }
  }
  return max_col+1;
}

/*! \brief Computes the degree of each feature, assuems degs init'd to all 0
 */
void CountDegrees(const vector::FVector<SVMExample> &ex, unsigned *degs) {
  for (size_t i = 0; i < ex.size; i++) {
    for (size_t j = 0; j < ex.values[i].vector.size; j++) {
      degs[ex.values[i].vector.index[j]]++;
    }
  }
}

}
}
}
