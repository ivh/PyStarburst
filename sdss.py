#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import pylab as P
import numpy as N
import PyGalKin.tool as T
import pyfits as F

def splitfits(fits):
    primHDU=fits[0]
    head=primHDU.header
    spec=primHDU.data[0,:]
    contsubspec=primHDU.data[1,:]
    noise=primHDU.data[2,:]
    mask=primHDU.data[3,:]
    sky=primHDU.data[4,:]

    return head,spec,noise

   
def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    
