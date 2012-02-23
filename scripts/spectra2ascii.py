#!/usr/bin/env python


import pyfits as F
import numpy as N
from sdss import *

OUTDIR='asciispec/'
conn,curs=DB.setupdb()

fields = 'objID,specObjID,plate,mjd,fiberID,z,Ha_w,Ha_h,Ha_s,Hb_w,Hd_w,O3727_s,agn,Mr,ext_g'

curs.execute('select %s from clean'%(fields,))

for values in curs:
    exec('%s = values'%(fields,))
    fname=DB.getspecfilename((mjd,plate,fiberID))
    try: fits=F.open(fname)
    except: stderr.write('Error in file %s\n'%fname); continue
    head,spec,noise=S.splitfits(fits)
    wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
    wave/=(1+z)
    outf=open(OUTDIR+'%s.spec'%specObjID,'w')
    outf.write('%d\n'%objID)
    outf.write('%s\n'%agn)
    outf.write( '\n'.join( ['%.8e'%x for x in Mr,ext_g,z,Ha_h,Ha_w,Ha_s,Hb_w,Hd_w,O3727_s] ) )
    outf.write('\n')
    for i,s in enumerate(spec):
            outf.write('%.8e %.8e %.8e\n'%(wave[i],s,noise[i]))
    outf.close()
