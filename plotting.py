#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import numpy as N
import pylab as P

from db import specfromid,fetch
from sdss import splitfits

def plotbyid(id):
    fits=specfromid(id)
    head,spec,noise=splitfits(fits)
    wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
    P.plot(wave,spec,linestyle='steps')


def plotHbs(cursor):
    cursor.execute('SELECT  2.5066*Hb_s*Hb_h/Hb_cont AS Hb_cw,Hb_w,p_15_w as Hb_pw FROM sdss WHERE zConf>0.95 AND Hb_w>0.1')
    Hb_cw,Hb_w,Hb_pw=fetch(cursor)
    P.subplot(221)
    P.plot(Hb_w,'.')
    P.title('Hb_w')
    P.subplot(222)
    P.plot(Hb_w/Hb_cw,'.')
    P.title('Hb_w/Hb_cw')
    P.subplot(223)
    P.plot(Hb_w/Hb_cw,'.')
    P.title('Hb_w/Hb_pw')
    P.subplot(224)
    P.plot(Hb_w/Hb_cw,'.')
    P.title('Hb_cw/Hb_pw')
    

def plotHas(cursor):
    cursor.execute('SELECT  2.5066*Ha_s*Ha_h/Ha_cont AS Ha_cw,Ha_w,p_15_w as Ha_pw FROM sb WHERE zConf>0.95 AND Ha_w>0.1 ORDER BY objID')
    Ha_cw,Ha_w,Ha_pw=fetch(cursor)
    P.subplot(221)
    P.plot(Ha_w,'.')
    P.title('Ha_w')
    P.subplot(222)
    P.plot(Ha_w/Ha_cw,'.')
    P.title('Ha_w/Ha_cw')
    P.subplot(223)
    P.plot(Ha_w/Ha_cw,'.')
    P.title('Ha_w/Ha_pw')
    P.subplot(224)
    P.plot(Ha_w/Ha_cw,'.')
    P.title('Ha_cw/Ha_pw')
    

   
def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    
