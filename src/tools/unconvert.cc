
// Offset to modify rows and columns by, may be set to 1 by compiler
// to subtract 1 from each row and column to compensate for matlab files
#ifndef OFFSET
#define OFFSET 0
#endif

#include <inttypes.h>

#include "hazy/types/tuple.h"
#include "hazy/scan/binfscan.h"

using hazy::scan::BinaryFileScanner; 

int main(int argc, char** argv) {
  if (argc != 3) {
    printf("usage: unconvert INFILE OUTFILE\n");
    printf("  converts binary to TSV, e.g. `convert in.bin out.tsv'\n");
    return 0;
  }
  assert(argc == 3);
  char *in = argv[1];
  char *out = argv[2];

  BinaryFileScanner scan(in);

  // now write it out
  FILE* f = fopen(out, "w");
  if (!f) {
    perror("cannot open output file");
    return 0;
  }

  while (scan.HasNext()) {
    hazy::types::Entry const &e = scan.Next();
    fprintf(f, "%d\t%d\t%lf\n", e.row+OFFSET, e.col+OFFSET, e.rating);
  }
  fclose(f);
}

