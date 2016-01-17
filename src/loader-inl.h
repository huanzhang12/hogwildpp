
#ifndef HAZY_LOADER_INL_H
#define HAZY_LOADER_INL_H

#include <vector>
#include <set>
#include <iostream>
#include <fstream>

#include "hazy/vector/fvector.h"

namespace hazy {
namespace hogwild {

std::vector<double>
load_floatlist(char *fname) {
  std::vector<double> s;
  std::ifstream f_in;
  double i;
  f_in.open(fname, std::ios::in);
  while(!f_in.eof()) {
    i = -1;
    f_in >> i;
    if (i == -1) {
      break;
    }
    s.push_back(i);
  }
  f_in.close();
  return s;
}


template <class Scan>
size_t LoadExamples(Scan &scan, vector::FVector<types::Entry> &fv) {
  scan.Reset();
  size_t max_col = 0;
  std::vector<types::Entry> vec;

  while (scan.HasNext()) {
    const types::Entry &e = scan.Next();
    if ((e.col > 0) && (static_cast<unsigned>(e.col) > max_col)) {
      max_col = e.col;
    }
    vec.push_back(e);
  }

  fv.values = new types::Entry[vec.size()];
  fv.size = vec.size();
  for (size_t i = 0; i < vec.size(); i++) {
    fv.values[i] = vec[i];
  }
  return max_col+1;
}

}
}
#endif
