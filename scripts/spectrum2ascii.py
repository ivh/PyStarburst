#!/usr/bin/env python


import pyfits as F
import numpy as N
from sys import stderr,argv,exit,stdout

def splitfits(fits):
    primHDU=fits[0]
    head=primHDU.header
    spec=primHDU.data[0,:]
    contsubspec=primHDU.data[1,:]
    noise=primHDU.data[2,:]
    mask=primHDU.data[3,:]
    sky=primHDU.data[4,:]

    return head,spec,noise

try:
    fname=argv[1]
    z=float(argv[2])
    fits=F.open(fname)
except:
    stderr.write('Error opening file or reading redshift.\n')
    exit()
head,spec,noise=splitfits(fits)
wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
wave/=(1+z)
for i,s in enumerate(spec):
    stdout.write('%.8e %.8e %.8e\n'%(wave[i],s,noise[i]))
  
