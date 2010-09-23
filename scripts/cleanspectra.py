#!/usr/bin/env python

USAGE="""Usage: ./getspectra.py <inputfile>
The input file must contain plate,mjd,fiberID."""

BASE='http://das.sdss.org/spectro/1d_26/'

import sys
from os import system, path, rename
from time import sleep

try: f=open(sys.argv[1])
except:
    print USAGE
    exit()

f.readline()

for line in f:
    specObjID,plate,mjd,fiberID,zConf,zWarning=line[:-1].split(',')
    fname='spSpec-%05d-%04d-%03d.fit'%(int(mjd),int(plate),int(fiberID))

    comm='mv toomany/%s spectra/'%fname
    system(comm)

