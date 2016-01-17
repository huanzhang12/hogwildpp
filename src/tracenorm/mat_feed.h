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


// Hogwild!, part of the Hazy Project
// research.cs.wisc.edu/hazy/
// Victor Bittorf, bittorf {at} cs.wisc.edu
// Chris Re, chrisre {at} cs.wisc.edu

#ifndef HAZY_HOGWILD_INSTANCES_MATFACT_MAT_FEED_H
#define HAZY_HOGWILD_INSTANCES_MATFACT_MAT_FEED_H
#include <cmath>
#include <iostream>
#include <cmath>

#include "ffa/instance.h"
#include "ffa/feedmodel.h"
#include "ffa/simplescan.h"
#include "ffa/memblockscan.h"
#include "hogwild/ffa/filescan.h"
#include "common/entrybuf.h"
#include "matfact/mat_model.h"
#include "matfact/mat_exec.h"

namespace hazy {
namespace hogwild {
namespace instances {
//! Tracenorm factorization
namespace matfact {


//! Creates a new Matrix Factorization Trough using the given memory
/*! Uses the provided buffer to allocate the trough and all memory needed
 * for double buffering and in-memory scanners. This could be a lot of
 * memory. 
 * \param mem large buffer used to back all allocations
 * \param fscan train file scanner
 * \param tscan test file scanner
 * \return A FeedTrough inialized and ready for ffa, backed by mem.
 */

template <class TrainFScan, class TestFScan>
class MFFeeder {
 public:
  typedef freeforall::FileScan<util::EntryBuffer, TrainFScan> FileScan_t;
  typedef freeforall::MemoryBlockScan<util::EntryBuffer> MemScan_t;
  typedef freeforall::Instance<Model_t, Params_t, MFExec> MFInstance;
  typedef freeforall::SimpleScan<util::EntryBuffer> TestScan_t;

  template <class TrainScan_t>
  static freeforall::FeedTrough<MFInstance, TrainScan_t, TestScan_t>* 
      Make(TrainFScan &fscan, TestFScan &tscan, Params_t &params, 
           size_t max_memory, const char *infile) {
    using namespace hazy;
    using mem::Buffer;
    using hogwild::util::EntryBuffer;
    size_t MByte = 1024*1024ULL;

    fscan.Reset();
    // pre-process each training example
    int nRows = 0;
    int nCols = 0;
    unsigned count = 0;
    while (fscan.HasNext()) {
      util::Entry const e = fscan.Next();
      if (nRows < e.row) nRows = e.row;
      if (nCols < e.col) nCols = e.col;
      count++;
    }

    size_t size_model = (nRows) * 4 * params.max_rank + (8+16) * nRows;
    size_model += (nCols) * 4 * params.max_rank + (8+16) * nCols;
    size_model += 1024;

    size_t size_params = nRows * 4 + nCols * 4;

    max_memory -= (size_params + size_model + MByte);

    if (max_memory <= 1024*1024*512) {
      printf("Not enough memory left! use -maxmem to give more memory\n");
      exit(-1);
    }

    fscan.Reset();
    printf("Rows = %d, Cols = %d\n", nRows, nCols);
    params.setup(nRows+1, nCols+1, count);

    while (fscan.HasNext()) {
      util::Entry const e = fscan.Next();
      parameter_map(params, e);
    }
    // Reset the scanner for the framework to use
    fscan.Reset();



    size_t block_size = max_memory / 2;
    printf("Page size is %.2fGB\n", block_size*1.0/(1024*1024*1024LLU));

    // Allocate our in-memory train scanner (double buffered)
    Buffer *exbuf_back = Buffer::Malloc(block_size);
    EntryBuffer *exbuf = new EntryBuffer(exbuf_back);
    Buffer *shadow_back = Buffer::Malloc(block_size);
    EntryBuffer *shadow = new EntryBuffer(shadow_back);

    // Make our test example buffer
    EntryBuffer *testbuf = new EntryBuffer(
        Buffer::Malloc(MByte*1024)); // XXX FIXME
    testbuf->Fill(tscan); // load testing examples, blocking call

    // Create the Scanners and start the initial permute
    TestScan_t *scan_test = new freeforall::SimpleScan<util::EntryBuffer>(testbuf);

    TrainScan_t *scan = new TrainScan_t(exbuf, shadow, fscan);
    scan->StartBasePermute();
    Benchmark::RunNCols::Value = params.nCols;
    Benchmark::RunNTrainEx::Value = count;
    Benchmark::RunNTestEx::Value = testbuf->GetExamples()->size;

    Model_t *m = new Model_t(params);
    MFInstance *inst = new MFInstance(m, &params);
    assert(params.L_degree != params.R_degree);

    if (infile != NULL) {
      m->LoadFromFile(infile);
    }

    return new freeforall::FeedTrough<MFInstance, TrainScan_t, TestScan_t>
        (inst, scan, scan_test);
  }

};

} // namespace matfact
} // namespace instances
} // namespace hogwild

} // namespace hazy

#endif


