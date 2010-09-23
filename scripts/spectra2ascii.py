#!/usr/bin/env python


import pyfits as F
import numpy as N
from sdss.db import getspecfilename
from sdss.sdss import splitfits
from sys import stderr,argv

datfile=open(argv[1])
datfile.readline()

for line in datfile:
    specObjID,plate,mjd,fiberID,z,Ha_w,Ha_h,Ha_s,Hb_w,Hd_w=line.split(',')
    plate,mjd,fiberID=map(int,(plate,mjd,fiberID))
    z,Ha_w,Ha_h,Ha_s,Hb_w,Hd_w=map(float,(z,Ha_w,Ha_h,Ha_s,Hb_w,Hd_w))

    fname=getspecfilename((mjd,plate,fiberID))
    try: fits=F.open(fname)
    except: stderr.write('Error in file %s\n'%fname); continue
    head,spec,noise=splitfits(fits)
    wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
    wave/=(1+z)
    outf=open('%s.spec'%specObjID,'w')
    outf.write('%.8e\n'%z)
    outf.write('%.8e\n'%Ha_h)
    outf.write('%.8e\n'%Ha_w)
    outf.write('%.8e\n'%Ha_s)
    outf.write('%.8e\n'%Hb_w)
    outf.write('%.8e\n'%Hd_w)
    for i,s in enumerate(spec):
            outf.write('%.8e %.8e %.8e\n'%(wave[i],s,noise[i]))
    outf.close()
  
