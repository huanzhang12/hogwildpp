// Copyright 2012 Victor Bittorf, Chris Re
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// The Hazy Project, http://research.cs.wisc.edu/hazy/
// Author : Victor Bittorf (bittorf [at] cs.wisc.edu)

#ifndef HAZY_CONCUR_BARRIER_T_H
#define HAZY_CONCUR_BARRIER_T_H

#include <pthread.h>

/* The prupose of this file is to allow the threading tools to be ported to
 * MacOSX which does not support pthread barriers. This is a naive barrier
 * implemtnation for the port which defaults to use pthreads when possible.
 */

namespace hazy {
namespace thread {
#ifdef __APPLE__
// we have to use our own barrier and timer
struct barrier_t {
  pthread_mutex_t mux;
  pthread_cond_t cond;
  int total;
  int current;
};
#else
typedef pthread_barrier_t barrier_t;
#endif



int barrier_init(barrier_t *b, void* attr, int count) {
#ifdef __APPLE__
  pthread_mutex_init(&b->mux, NULL);
  pthread_cond_init(&b->cond, NULL);
  b->total = count;
  b->current = count;
  return 0;
#else
  return pthread_barrier_init(b, (pthread_barrierattr_t*)attr, count);
#endif
}

int barrier_wait(barrier_t *b) {
#ifdef __APPLE__
  pthread_mutex_lock(&b->mux);
  b->current--;
  if (b->current == 0) {
    // reset the count
    b->current = b->total;
    // wake up everyone, we're the last to the fence
    pthread_cond_broadcast(&b->cond);
    pthread_mutex_unlock(&b->mux);
    return 0;
  } 
  // otherwise, wait the fence
  pthread_cond_wait(&b->cond, &b->mux);
  // release the mux
  pthread_mutex_unlock(&b->mux);
  return 0;
#else
  return pthread_barrier_wait(b);
#endif
}


int barrier_destroy(barrier_t *b) {
#ifdef __APPLE__
  // XXX FIXME TODO
  return -1;
#else
  return pthread_barrier_destroy(b);
#endif
}

} // namespace thread
} // namespace hazy
#endif
