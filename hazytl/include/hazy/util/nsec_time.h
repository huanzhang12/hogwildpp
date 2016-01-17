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
//
/* This file gives access to the current nano-second time
 */
#ifndef HAZY_UTIL_NSEC_TIME_H
#define HAZY_UTIL_NSEC_TIME_H

#include <time.h>
#include <sys/time.h>

#ifdef APPLE
#include <mach/clock.h>
#include <mach/mach.h>
#endif

namespace hazy {
namespace util {

/*! \brief Get the current nano-seconds since some epoch
 * This is portable to OSX as well as Linux
 * \return nsec since epoch
 */
inline unsigned long long CurrentNSec() {
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec * 1000 *1000 *1000UL + tv.tv_usec*1000UL;
#else
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return (unsigned long long) ts.tv_nsec + 
      (unsigned long long) ts.tv_sec * 1000 * 1000 * 1000;
#endif
}

} // namespace util
} // namespace hazy
#endif
