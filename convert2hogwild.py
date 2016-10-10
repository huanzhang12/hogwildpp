#!/usr/bin/env pypy

import os, sys

if len(sys.argv) < 2:
	print "Usage: {} <input binary> <output tsv>".format(sys.argv[0])
	os._exit(0)

input_file = sys.argv[1]
output_file = sys.argv[2]

out = open(output_file, 'w')
index = 0

print "Processing {}".format(input_file)
with open(input_file, 'r') as f:
	for line in f:
		elements = line.split()
		print >> out, "{}\t{}\t{}".format(index, -2, elements[0])
		for values in elements[1:]:
			pairs = values.split(':')
			print >> out, "{}\t{}\t{}".format(index, int(pairs[0])-1, pairs[1])
		index += 1
		if (index & 0xfff == 0xfff):
			sys.stdout.write('.')
			sys.stdout.flush()


print '\n{} lines processed'.format(index)

