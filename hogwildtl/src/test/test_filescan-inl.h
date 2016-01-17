#include "gtest/gtest.h"

#include "hazy/scan/tsvfscan.h"
#include "hazy/types/tuple.h"

#include "hazy/hogwild/hogwild-inl.h"
#include "hazy/hogwild/file_scan.h"

using namespace hazy;

TEST(TSVScan, SimpleLoad) {
  scan::TSVFileScanner scan("test_data/simple.tsv");

  ASSERT_TRUE(scan.HasNext());

  unsigned r = 0;
  while (scan.HasNext()) {
    const types::Entry &ent  = scan.Next();
    ASSERT_LT(r, 8u);
    ASSERT_EQ((int)(r >> 2), ent.row);
    ASSERT_EQ((int)(r & 0x3), ent.col);
    ASSERT_EQ(r, ent.rating);
    r++;
  }
  ASSERT_EQ(r, 8u);
}

TEST(TSVScan, TwoPasses) {
  scan::TSVFileScanner scan("test_data/simple.tsv");

  ASSERT_TRUE(scan.HasNext());

  unsigned r = 0;
  while (scan.HasNext()) {
    const types::Entry &ent  = scan.Next();
    ASSERT_LT(r, 8u);
    ASSERT_EQ((int)(r >> 2), ent.row);
    ASSERT_EQ((int)(r & 0x3), ent.col);
    ASSERT_EQ(r, ent.rating);
    r++;
  }
  ASSERT_FALSE(scan.HasNext());

  scan.Reset();

  ASSERT_TRUE(scan.HasNext());

  r = 0;
  while (scan.HasNext()) {
    const types::Entry &ent  = scan.Next();
    ASSERT_LT(r, 8u);
    ASSERT_EQ((int)(r >> 2), ent.row);
    ASSERT_EQ((int)(r & 0x3), ent.col);
    ASSERT_EQ(r, ent.rating);
    r++;
  }
  ASSERT_EQ(r, 8u);

}




TEST(TSVScan, MidReset) {
  scan::TSVFileScanner scan("test_data/simple.tsv");

  ASSERT_TRUE(scan.HasNext());
  scan.Next();

  scan.Reset();

  unsigned r = 0;
  while (scan.HasNext()) {
    const types::Entry &ent  = scan.Next();
    ASSERT_LT(r, 8u);
    ASSERT_EQ((int)(r >> 2), ent.row);
    ASSERT_EQ((int)(r & 0x3), ent.col);
    ASSERT_EQ(r, ent.rating);
    r++;
  }
  ASSERT_EQ(r, 8u);
}



TEST(FileScan, SmallBuffer) {
  scan::TSVFileScanner scan("test_data/simple.tsv");

  hogwild::FileScan<scan::TSVFileScanner, types::Entry> fscan(scan, 64);
  fscan.Init();

  ASSERT_TRUE(fscan.HasNext());

  unsigned r = 0;
  while (fscan.HasNext()) {
    hogwild::ExampleBlock<types::Entry> &blk = fscan.Next();
    for (unsigned i = 0; i < blk.ex.size; i++) {
      ASSERT_LT(r, 8u);
      ASSERT_EQ((int)(r >> 2), blk.ex.values[i].row);
      ASSERT_EQ((int)(r & 0x3), blk.ex.values[i].col);
      ASSERT_EQ(r, blk.ex.values[i].rating);
      r++;
    }
  }
  ASSERT_EQ(r, 8u);

  fscan.Destroy();
}


TEST(FileScan, TwoScans) {
  scan::TSVFileScanner scan("test_data/simple.tsv");

  hogwild::FileScan<scan::TSVFileScanner, types::Entry> fscan(scan, 32);
  fscan.Init();

  ASSERT_TRUE(fscan.HasNext());

  unsigned r = 0;
  while (fscan.HasNext()) {
    hogwild::ExampleBlock<types::Entry> &blk = fscan.Next();
    for (unsigned i = 0; i < blk.ex.size; i++) {
      ASSERT_LT(r, 8u);
      ASSERT_EQ((int)(r >> 2), blk.ex.values[i].row);
      ASSERT_EQ((int)(r & 0x3), blk.ex.values[i].col);
      ASSERT_EQ(r, blk.ex.values[i].rating);
      r++;
    }
  }
  ASSERT_EQ(r, 8u);

  fscan.Reset();
  r = 0;
  while (fscan.HasNext()) {
    hogwild::ExampleBlock<types::Entry> &blk = fscan.Next();
    for (unsigned i = 0; i < blk.ex.size; i++) {
      ASSERT_LT(r, 8u);
      ASSERT_EQ((int)(r >> 2), blk.ex.values[i].row);
      ASSERT_EQ((int)(r & 0x3), blk.ex.values[i].col);
      ASSERT_EQ(r, blk.ex.values[i].rating);
      r++;
    }
  }
  ASSERT_EQ(r, 8u);
  fscan.Destroy();
}

TEST(FileScan, BigBuffer) {
  scan::TSVFileScanner scan("test_data/simple.tsv");

  hogwild::FileScan<scan::TSVFileScanner, types::Entry> fscan(scan, 1024);
  fscan.Init();

  ASSERT_TRUE(fscan.HasNext());

  unsigned r = 0;
  hogwild::ExampleBlock<types::Entry> &blk = fscan.Next();
  for (unsigned i = 0; i < blk.ex.size; i++) {
    ASSERT_LT(r, 8u);
    ASSERT_EQ((int)(r >> 2), blk.ex.values[i].row);
    ASSERT_EQ((int)(r & 0x3), blk.ex.values[i].col);
    ASSERT_EQ(r, blk.ex.values[i].rating);
    r++;
  }
  ASSERT_EQ(r, 8u);

  fscan.Destroy();
}

