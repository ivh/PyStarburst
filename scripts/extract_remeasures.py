#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
    print 'give input files as arguments!'
    sys.exit()
else:
    infiles=sys.argv[1:]

outf = open('remeas_ew','w')
outf.write('SpecID Ha Hg Hd He\n')

for fname in infiles:
    outf.write(fname.split('/')[-1])
    for line in open(fname):
        if line.startswith('SDSS EW'):
            val = line.split(':')[1].strip()
            outf.write(' '+val)
        elif line.startswith('Model'):
            break
    outf.write('\n')
