#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import pylab as P
import numpy as N
import pyfits as F

# CONSTANTS
#
# units:
#        velocity: km/s
#        wavelength: Angstrom
lambHA=N.array([6562.7797852000003],'Float32')
sol=N.array([299792.458],'Float32')
c=sol
H0=N.array([72.],'Float32')
Grav=N.array([6.6726E-11*1.989E30/1000**3],'Float32') ### in solar masses and km
pc=N.array([3.086E13],'Float32') ## in km


def splitfits(fits):
    primHDU=fits[0]
    head=primHDU.header
    spec=primHDU.data[0,:]
    contsubspec=primHDU.data[1,:]
    noise=primHDU.data[2,:]
    mask=primHDU.data[3,:]
    sky=primHDU.data[4,:]

    return head,spec,noise


def calcHaW(Ha_cont,Ha_h,Ha_s):
    return 2.5066 * Ha_s * Ha_h / Ha_cont

def Ang2KMS(ang,z,restlamb=lambHA):
    return ang/(1+z)/restlamb*sol

def fixSigmaInDB(curs,which='Ha_s',restlamb=lambHA):
    query=curs.execute("SELECT objID,%s,z from sdss WHERE (%s NOTNULL) AND (z NOTNULL) "%(which,which))
    all=query.fetchall()
    for id,sig,z in all:
        curs.execute("UPDATE sdss SET %s=%e WHERE (objID=%s)"%(which,Ang2KMS(sig,z,restlamb),id))
    print "Dont forget to commit!"
    
def extcorr(ec,filt):
    if filt == 'u': const=0.7156
    elif filt == 'g': const=0.5031
    elif filt == 'r': const=0.3518
    elif filt == 'i': const=0.3013 
    elif filt == 'z': const=0.2313
    else:
        print "Unknown filter"
        return False
    return 2.5*const*ec

def fillExtCorrDB(curs):
    filters=['u','g','r','i','z']
    for filt in filters:
        
        print extcorr

def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    
