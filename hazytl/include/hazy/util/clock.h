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

#ifndef HAZY_UTIL_CLOCK_H
#define HAZY_UTIL_CLOCK_H

#include "hazy/util/nsec_time.h"

namespace hazy {
namespace util {

/*! \brief A simple timer which can be started and stopped
 * Only read Clock.value after a stop or pause, otherwise use Read()
 */
class Clock {
 public:
  double value; //< the elapsed time in seconds

  Clock() : value(0), start(0), stop(0),  paused_(false) { } 

  /*! \brief Start counting.
   */
  inline void Start() {
    start = CurrentNSec();
    if (!paused_) {
      value = 0.0;
    }
  }

  /*! \brief Return the delta from the last call to Start()
   */
  inline double Read() {
    if (paused_) {
      return value;
    } else {
      stop = CurrentNSec();
      value = (double) (stop - start) * 1e-9;
      return value;
    }
  }

  /*! \brief Pause counting, keep a running sum of time
   */
  inline double Pause() {
    stop = CurrentNSec();
    value += (double) (stop - start) * 1e-9;
    paused_ = true;
    return value;
  }

  /*! \brief Stop counting, next call to Start() will reset the time
   */
  inline double Stop() {
    stop = CurrentNSec();
    value = (double) (stop - start) * 1e-9;
    paused_ = false;
    return value;
  }

 private:
  unsigned long long start;
  unsigned long long stop;
  bool paused_;
};

} // namespace util
} // namespace hazy
#endif
