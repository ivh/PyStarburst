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
from sdss import *

def micJy2Watt(mJy,z,lambd):
    lambd*=1E-10
    dis=sdss.distanceInMeter(z)
    a=4*pi*(dis**2)
    return mJy*1E-32*a*(c*1000)/lambd

def micJy2SolarLum(mJy,z,lambd=1530):
    solarLum=3.846E26 #Watt
    return micJy2Watt(mJy,z,lambd)/solarLum

def uext2fuv(uext):
    """ gives the factor to multiply the fuv-flux with, dependign on u-band extiction in magnitudes"""
    return 1.8 * (10**(uext/2.5))

def selectdata():
    """retrun the query"""
    query="SELECT gal.objID, spec.specobjID, gal.ra, gal.dec, spec.z"
    query+=", gal.extinction_u ext_u"
    query+=", Ha.ew Ha_w, Ha.height Ha_h, Ha.sigma Ha_s"
    query+=", Hb.height Hb_h"
    query+=", OIII.height OIII_h"
    query+=", NII.height NII_h"
    

    query+=" INTO mydb.galex2"
    query+=" FROM mydb.galex1 sel, Galaxy gal, SpecObj spec, SpecLine Ha, SpecLine Hb, SpecLine OIII, SpecLine NII"

    query+=" WHERE (gal.objID=sel.objID)"
    query+=" AND (spec.bestObjID=sel.objID)"
    query+=" AND (Ha.specobjID=spec.specobjID)"
    query+=" AND (Hb.specobjID=spec.specobjID)"
    query+=" AND (OIII.specobjID=spec.specobjID)"
    query+=" AND (NII.specobjID=spec.specobjID)"
    
    query+=" AND (Ha.LineID=6565)"
    query+=" AND (Hb.LineID=4863)" 
    query+=" AND (OIII.LineID=5008)" 
    query+=" AND (NII.LineID=6585)"

    return query

def get_sdss(file='/home/tom/galex/galex2_thomasmarquart.csv'):
    objID,specobjID=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Int64',usecols=[0,1]))
    ra,dec,z,ext_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Float64',usecols=[2,3,4,5,6,7,8,9,10,11]))
    return objID,specobjID,ra,dec,z,ext_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h

def get_galex(file='/home/tom/galex/fluxes_thomasmarquart.csv'):
    sid,gid=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Int64',usecols=[0,1]))
    ra,dec,z,fuv_flux,nuv_flux,fuv_corr=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Float64',usecols=[2,3,4,5,6,7]))
    return sid,gid,ra,dec,z,fuv_flux,nuv_flux,fuv_corr
    
def plotAGN():
    objID,specobjID,ra,dec,z,ext_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h=get_sdss()
    
    x=P.log10(NII_h/Ha_h)
    y=P.log10(OIII_h/Hb_h)
    xl=P.arange(-1.7,-0.15,0.01)
    yl=sdss.mylee(xl)
    P.plot(x,y,',',label='all 2011 objects')
    broadx=N.where(Ha_s>5,N.log10(NII_h/Ha_h),0.0)
    broady=N.where(Ha_s>5,N.log10(OIII_h/Hb_h),0.0)
    P.plot(broadx,broady,'r,',label=r'sigma(H-alpha) > 5$\AA$')

    P.plot(xl,yl,'r-',linewidth=2)
    P.legend()
    P.xlabel('log [NII]/Ha')
    P.ylabel('log [OIII]/Hb')

def correct_ext():
    conn,curs=setupdb('/home/tom/galex/test.db')
    query=curs.execute("SELECT g.sid,g.gid,s.ext_u,g.fuv_flux FROM sdss s, galex g WHERE g.sid=s.objID")
    for sid,gid,ext,fuv in query.fetchall():
        que="UPDATE galex SET fuv_corr=%f WHERE gid=%d"%(fuv*uext2fuv(ext),gid)
        curs.execute(que)
    conn.commit()
    conn.close()

def plot_z_fuv():
    conn,curs=setupdb('/home/tom/galex/test.db')
    z,fuv=gettable(curs,'s.z,g.fuv_corr',table='sdss s, galex g',where='g.sid=s.objID AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=1')
    P.semilogy(z,micJy2SolarLum(fuv,z),'r,')
    P.xlabel('z')
    P.ylabel(r'L$_{FUV}$')
    conn.close()

def decideAGN():
    conn,curs=setupdb('/home/tom/galex/test.db')
    ids=gettable(curs,cols='objid',where='(Ha_h >0) AND (Hb_h>0) AND (OIII_h>0) AND (NII_h>0)',table='sdss')[0]
    x,y,sig=gettable(curs,cols='NII_h/Ha_h,OIII_h/Hb_h,Ha_s',where='(Ha_h >0) AND (Hb_h>0) AND (OIII_h>0) AND (NII_h>0)',table='sdss')
    x=N.log10(N.array(x))
    y=N.log10(N.array(y))
    limit=mylee(x)
    for i,id in enumerate(ids):
        #print id,y[i],limit[i]
        if (y[i]<limit[i]) and (sig[i] < 5):
            curs.execute("UPDATE sdss SET agn=0 WHERE objid=%s"%(id,))
        else:
            print x[i],y[i],limit[i],sig[i]
            curs.execute("UPDATE sdss SET agn=1 WHERE objid=%s"%(id,))

    conn.commit()
    conn.close()