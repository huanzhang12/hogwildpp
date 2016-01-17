
#ifndef HAZY_HOGWILD_HOGWILD_TASK_H
#define HAZY_HOGWILD_HOGWILD_TASK_H

#include "hazy/vector/fvector.h"

namespace hazy {
namespace hogwild {

template <class Example>
struct ExampleBlock {
    vector::FVector<Example> ex;
    vector::FVector<size_t> perm;
};

template <class Model, class Params, class Example>
struct HogwildTask {
  Model *model;
  Params *params;
  ExampleBlock<Example> *block;
};

} // namespace hogwild
} // namespace hazy
#endif
