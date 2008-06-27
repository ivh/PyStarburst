#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import numpy as N
import pylab as P

from db import specfromid

def plotbyid(id):
    fits=specfromid(id)
    P.plot(fits[0].data[0,:])


   
def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    
