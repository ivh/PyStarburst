#!/usr/bin/env python

USAGE="""Usage: ./getspectra.py <inputfile>
The input file must contain plate,mjd,fiberID."""

BASE='http://das.sdss.org/spectro/1d_26/'

import sys
from os import system, path
from time import sleep

try: f=open(sys.argv[1])
except:
    print USAGE
    exit()

#f.readline()

for line in f:
    specObjID,plate,mjd,fiberID=line.strip().split(',')
    fname='spSpec-%05d-%04d-%03d.fit'%(int(mjd),int(plate),int(fiberID))

    if path.exists(fname): print '%s already here'%fname; continue

    print specObjID,plate,mjd,fiberID,' ... ',
    comm='wget -c -q %s%04d/1d/%s'%(BASE,int(plate),fname)
    system(comm)
    print 'done'
    sleep(0.5)
