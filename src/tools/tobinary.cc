
#include <inttypes.h>

#include "hazy/scan/tsvfscan.h"
#include "hazy/types/tuple.h"

// Offset to modify rows and columns by, may be set to 1 by compiler
// to subtract 1 from each row and column to compensate for matlab files
#ifndef MATLAB_CONVERT_OFFSET
#define MATLAB_CONVERT_OFFSET 0
#endif

using hazy::scan::TSVFileScanner; 

int main(int argc, char** argv) {
  if (argc != 3) {
    printf("usage: convert INFILE OUTFILE\n");
    printf("  converts TSV to binary, e.g. `convert in.tsv out.bin'\n");
    return 0;
  }
  assert(argc == 3);
  char *in = argv[1];
  char *out = argv[2];

  TSVFileScanner scan(in);

  // first count how many there are
  uint64_t count = 0;
  while (scan.HasNext()) {
    scan.Next();
    count++;
  }

  // now write it out
  FILE* f = fopen(out, "w");
  if (!f) {
    perror("cannot open output file");
    return 0;
  }

  printf("Found %lu examples.\n", count);
  fwrite(&count, sizeof(count), 1, f);



  hazy::types::Entry offset;
  TSVFileScanner scan2(in);
  while (scan2.HasNext()) {
    hazy::types::Entry const &e = scan2.Next();
    offset.rating = e.rating;
    offset.row = e.row - MATLAB_CONVERT_OFFSET;
    offset.col = e.col - MATLAB_CONVERT_OFFSET;

    fwrite(&offset, sizeof(offset), 1, f);
  }

  fclose(f);
}
