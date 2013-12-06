#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import pylab as P
import numpy as N
import pyfits as F
from db import *
from numpy.ma import masked_where
from math import pi

# CONSTANTS
#
# units:
#        velocity: km/s
#        wavelength: Angstrom
lambHA=N.array([6562.7797852000003],'Float32')
sol=N.array([299792.458],'Float32')
c=sol
H0=N.array([70.],'Float32')
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

def Kms2Ang(kms,z,restlamb=lambHA):
    return kms*(1+z)*restlamb/sol

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

def absmag(m,z):
    return m-(5*N.log10((z*3E5)/72))-25

def stasinska(x):
    return (-30.787 + (1.1358 * x) + (0.27297 * (x**2))) * N.tanh(5.7409 * x) - 31.093

def kauffman(x):
    return 0.61/(x-0.05)+1.3

def lee(x):
    return 0.61/(x-0.25)+1.25

def mylee(x):
    return 0.71/(x-0.25)+1.25

def schechterBlanton(m):
    """
    returns the luminosity function on unit magnitudes
    as blanton etal http://iopscience.iop.org/1538-3881/121/5/2358/fulltext
    fitted it to sdss data
    """
    phi0 = 1.46E-2 # +- 0.12E-2; units h^3 / Mpc^3
    phi0 *= 0.7**3
    M0 = -20.83 # +. 0.03
    alpha = -1.20 # +- 0.03
    phi = 0.4 * N.log(10) * phi0
    phi *= 10**(-0.4*(m-M0)*(alpha+1))
    phi *= N.exp(-10**(-0.4*(m-M0)))

    return phi

def decideAGN(curs):
    ids,x,y,sig=gettable(curs,cols='objid,NII_h/Ha_h,O5008_h/Hb_h,Ha_s',where='(Ha_h >0) AND (Hb_h>0) AND (O5008_h>0) AND (NII_h>0)',table='sdss')
    x=N.log10(x.astype('float64'))
    y=N.log10(y.astype('float64'))
    limit=mylee(x)
    createcolumnifnotexists(curs,'agn')
    for i,id in enumerate(ids):
        #print id,y[i],limit[i]
        if (y[i]>limit[i]) or (sig[i] > 250) or (x[i]>0):
            curs.execute("UPDATE sdss SET agn=1 WHERE objid=%s"%(id,))
        else:
            #print x[i],y[i],limit[i],sig[i]
            curs.execute("UPDATE sdss SET agn=0 WHERE objid=%s"%(id,))

def fillExtCorrDB(curs):
    filters=['u','g','r','i','z']
    ec=gettable(curs,'ec','ec NOTNULL')
    for filt in filters:

        print extcorr(ec,filt)


def fillM(curs):
    createcolumnifnotexists(curs,'Mr')
    createcolumnifnotexists(curs,'Mg')
    ids,m,z=gettable(curs,cols='objid,m_g-ext_g,z',where='m_g NOTNULL AND z NOTNULL',table='sdss')
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET Mg=%f WHERE objid=%s"%(absmag(m[i],z[i]),id))
    ids,m,z=gettable(curs,cols='objid,m_r-ext_r,z',where='m_r NOTNULL AND z NOTNULL',table='sdss')
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET Mr=%f WHERE objid=%s"%(absmag(m[i],z[i]),id))

def massG(m,ml):
    return ml*10**(-0.4*(m-5.48))

def fillMass(curs,table='sbfit',ml='ML'):
    createcolumnifnotexists(curs,'mass',table=table)
    curs.execute('select f.ID,f.%s,s.f_g,s.ext_g,s.z from sdss s, %s f where s.specObjID=f.ID'%(ml,table))
    data=curs.fetchall()
    for id,ml,mg,ext,z in data:
        curs.execute("UPDATE %s SET mass=%f WHERE ID=%s"%(table,massG(absmag(mg-ext,z),ml),id))

def z2volume(z):
    return (4*P.pi/3)*(z*c/72.)**3

def voldens(Mr):
    r=10**(0.2*(17.77-Mr+5))
    r/=1E6
    v=4.0 * P.pi * (r**3) / 3.0
    v/=5.13 # SDSS DR7 spectra cover 19% of the sky
    v-=z2volume(0.005) # subtract the local volume
    return 1/v

def fillVoldens(curs):
    createcolumnifnotexists(curs,'voldens')
    ids,r,z=gettable(curs,cols='objid,m_r,z',where='m_r NOTNULL AND z NOTNULL',table='sdss')
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET voldens=%e WHERE objid=%s"%(voldens(absmag(r[i],z[i])),id))

def fadingR(age):
    if age<=1E7: return 3.6
    if age>1E7 and age<=1E9: return 3.6-(1.6*(N.log10(age)-7))
    else: return 0.0

def fillFading(curs):
    createcolumnifnotexists(curs,'fade',table='sbfit')
    ids,age,z=gettable(curs,cols='ID,Age,z',where='age NOTNULL',table='sbfit')
    for i,id in enumerate(ids):
        curs.execute("UPDATE sbfit SET fade=%f WHERE ID=%s"%(fadingR(age[i]),id))

def distanceInMeter(z):
    dist=(z * c / H0) * 3.0856776e+22
    return dist.astype('Float64')

def arcsec2meter(arcsec,z):
    D=distanceInMeter(z)
    a=N.radians(arcsec/3600.)
    return a*D

def arcsec2kpc(arcsec,z):
    return arcsec2meter(arcsec,z) / 3.0856776e19

def sfr(Ha_h,Ha_s,z,ec=0,corrFac=1):
    #            ||  extra -3 because of w/m^2 vs erg/s/cm^2
    return corrFac*Ha_h*1E-20 *Kms2Ang(Ha_s,z)*N.sqrt(2*P.pi) *4*P.pi*distanceInMeter(z)**2 / 1.51e34 * (10**(0.4*extcorr(ec,'r')))

def fillSFR(curs):
    createcolumnifnotexists(curs,'sfr',table='sbfit')
    data=curs.execute('select f.ID,f.log_Haflux,s.z,s.Ha_h,s.Ha_s from sdss s, sbfit f where s.specObjID=f.ID and s.Ha_h >0 AND s.Ha_s>1 AND f.log_Haflux NOTNULL').fetchall()
    for id,Ha_lum,z,Ha_h,Ha_s in data:
        curs.execute("UPDATE sbfit SET sfr=%f WHERE ID=%s"%(10**Ha_lum / 1.51e34,id))

def mgas(Mg):
    return 10**(-0.3536*Mg+2.6374) + 10**(-0.45*Mg+0.35)

def fillMgas(curs):
    createcolumnifnotexists(curs,'mgas')
    ids,m_g,z=gettable(curs,cols='objid,f_g,z',where='m_g NOTNULL AND z NOTNULL',table='sdss')
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET mgas=%f WHERE objid=%s"%(mgas(absmag(m_g[i],z[i])),id))

def fillMtot(curs,table='sbfit'):
    createcolumnifnotexists(curs,'mtot',table=table)
    data=curs.execute('select f.id,f.mass,s.mgas from %s f, sdss s where f.ID=s.specObjID'%table).fetchall()
    for id,mass,mgas in data:
        curs.execute("UPDATE %s SET mtot=%f WHERE ID=%s"%(table,mgas+mass,id))

def fillBpara2(curs):
    createcolumnifnotexists(curs,'bpara2',table='sbfit')
    ids,mass,sfr=gettable(curs,cols='ID,mass,sfr',where='mass NOTNULL AND sfr NOTNULL',table='sbfit')
    for i,id in enumerate(ids):
        curs.execute("UPDATE sbfit SET bpara2=%f WHERE ID=%s"%(sfr[i]/(mass[i]/1E10),id))

def lumfu(X,M,voldens):
    y=N.zeros(len(X),dtype='Float64')
    for i,x in enumerate(X):
        y[i] = N.sum(N.where(M<x,voldens,0.0))
    y[1:] = y[1:] - y[:-1]
    return y

def averbins(X,orgX,Y,median=True):
    mean=N.zeros(len(X),dtype='Float64')
    #median=mean.copy()
    #sigma=mean.copy()
    for i,x in enumerate(X):
        tmp=masked_where(orgX<x,Y)
        if i<len(X)-1: tmp=masked_where(orgX>=X[i+1],tmp)
        print i,x,tmp
        if median: mean[i]=N.ma.median(tmp)
        else: mean[i]=tmp.mean()
        #median[i]=N.median(tmp)
        #sigma[i]=tmp.std()
    return mean

def corr_sigma_instr(sigma,instrument_sigma=70):
    result = N.sqrt(sigma**2 - instrument_sigma**2)
    return N.where(N.isnan(result), 0.0, result)

def dynMassDisk(r,sigma):
    'r in kpc, sigma in km/s, returns solar masses'
    return 7.9E5 * r * corr_sigma_instr(sigma)**2

def dynMassSphere(r,sigma):
    'r in kpc, sigma in km/s, returns solar masses'
    return 1.1E6 * r * corr_sigma_instr(sigma)**2

def fillDynMasses(curs):
    createcolumnifnotexists(curs,'dynMassDisk')
    createcolumnifnotexists(curs,'dynMassSphere')
    ids,r50,r90,z,sig = gettable(curs,cols='objid,petroR50_r,petroR90_r,z,Ha_s',where='z NOTNULL AND Ha_s NOTNULL',table='sdss')
    ms = dynMassSphere(arcsec2kpc(r50,z),sig)
    md = dynMassDisk(arcsec2kpc(r90,z),sig)
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET dynMassDisk=%s, dynMassSphere=%s WHERE objid=%s"%(md[i],ms[i],id))

def fillBpara3(curs):
    createcolumnifnotexists(curs,'bpara3',table='sbfit')
    data=curs.execute('select f.ID,f.sfr,s.dynMassDisk,s.dynMassSphere from sbfit f, sdss s where f.ID=s.specObjID and f.sfr NOTNULL').fetchall()
    for id,sfr,m1,m2 in data:
        m=(m1+m2)/2.0
        try: b = sfr/(m/1E10)
        except: b='NULL'
        curs.execute("UPDATE sbfit SET bpara3=%s WHERE ID=%s"%(b,id))

def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()

