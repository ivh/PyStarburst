#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import numpy as N
import pylab as P

from db import *
from sdss import splitfits

def plotbyid(id):
    fits=specfromid(id)
    head,spec,noise=splitfits(fits)
    wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
    P.plot(wave,spec,linestyle='steps')


def plotHas(cursor):
    Ha_cw,Ha_w,Ha_pw,Ha_s=get(cursor,'SELECT  2.5066*Ha_s*Ha_h/Ha_cont AS Ha_cw,Ha_w,p_15_w as Ha_pw,Ha_s FROM sb WHERE zConf>0.95 AND Ha_w>0.1 AND Ha_w!=0.0 AND Ha_s BETWEEN 5 AND 15 ORDER BY Ha_w')
    P.subplot(221)
    P.plot(Ha_w,Ha_cw,'.')
    P.title('Ha_w vs. Ha_cw')
    P.subplot(222)
    P.plot(Ha_w,Ha_cw/Ha_w,'.')
    P.title('Ha_cw/Ha_w')
    P.subplot(223)
    P.plot(Ha_w,Ha_pw/Ha_w,'.')
    P.title('Ha_pw/Ha_w')
    P.subplot(224)
    P.plot(Ha_s,Ha_cw/Ha_w,'.')
    P.title('Ha_cw/Ha_w vs. sigma')
    

def plot1(cursor):
    big=getsb(cursor,'Ha_s,Ha_h,NII_h','Ha_s >15')
    small=getsb(cursor,'Ha_s,Ha_h,NII_h','Ha_s <5')
    midd=getsb(cursor,'Ha_s,Ha_h,NII_h','(Ha_s >5 AND Ha_s <15)')
    numsb=midd[0].size+big[0].size+small[0].size
    P.plot(N.log10(big[2]/big[1]),N.log10(big[0]),',k')
    P.plot(N.log10(small[2]/small[1]),N.log10(small[0]),',b')
    P.plot(N.log10(midd[2]/midd[1]),N.log10(midd[0]),',g')
    P.title('Starbursts (%d)'%numsb)
    P.xlabel('log NII/Ha')
    P.ylabel('log sigma Ha')

def plot2(cursor):
    big=getsb(cursor,'O5008_h/Hb_h,Ha_h,NII_h','Ha_s >15')
    small=getsb(cursor,'O5008_h/Hb_h,Ha_h,NII_h','Ha_s <5')
    midd=getsb(cursor,'O5008_h/Hb_h,Ha_h,NII_h','(Ha_s >5 AND Ha_s <15)')
    numsb=midd[0].size+big[0].size+small[0].size
    P.plot(N.log10(big[2]/big[1]),N.log10(big[0]),',k')
    P.plot(N.log10(small[2]/small[1]),N.log10(small[0]),',b')
    P.plot(N.log10(midd[2]/midd[1]),N.log10(midd[0]),',g')
    P.title('Starbursts, %s galaxies'%numsb)
    P.xlabel('log NII/Ha')
    P.ylabel('log O5008 / Hb')

    
def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    
