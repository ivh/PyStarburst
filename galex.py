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
from numpy import pi
from sdss import *



def micJy2Watt(mJy,z,lambd):
    lambd*=1E-10
    dis=sdss.distanceInMeter(z)
    a=4*pi*(dis**2)
    return mJy*1E-32*a*(c*1000)/lambd

def sdssflux2Watt(height,width,z):
    """ fluxes (the amplitude of the fitted gaussian) in sdss come in 1E-17 erg/s/cm**2/Å """
    dis=sdss.distanceInMeter(z)*100 # cm
    area=4*pi*(dis**2)
    return N.sqrt(2*pi)*1E-24*area*height*width # erg is 1E-7 J

def micJy2SolarLum(mJy,z,lambd=1530):
    solarLum=3.846E26 #Watt
    return micJy2Watt(mJy,z,lambd)/solarLum

def uext2fuv(uext):
    """ gives the factor to multiply the fuv-flux with, dependign on u-band extiction in magnitudes"""
    return 10**(1.66*uext/2.5)


#
# PLOTTING
#
   
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

def plot_z_fuv():
    conn,curs=setupdb('/home/tom/galex/data.db')
    z,fuv=gettable(curs,'s.z,g.fuv_lum',table='sdss s, galex g',where='g.sid=s.objID AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0')
    P.semilogy(z,fuv,'b,')
    P.xlabel('z')
    P.ylabel(r'L$_{FUV}$')
    conn.close()


#
# FILLING CERTAIN DATABASE-TABLES
#

def fill_ext():
    conn,curs=setupdb('/home/tom/galex/data.db')
    query=curs.execute("SELECT g.sid,g.gid,s.extinction_u,g.fuv_flux FROM sdss s, galex g WHERE g.sid=s.objID")
    for sid,gid,ext,fuv in query.fetchall():
        que="UPDATE galex SET fuv_corr=%f WHERE gid=%d"%(fuv*uext2fuv(ext),gid)
        curs.execute(que)
    conn.commit()
    conn.close()

def fill_fuv_lum():
    conn,curs=setupdb('/home/tom/galex/data.db')
    query=curs.execute("SELECT g.sid,g.gid,s.extinction_u,g.fuv_corr,s.z FROM sdss s, galex g WHERE g.sid=s.objID")
    for sid,gid,ext,fuv,z in query.fetchall():
        que="UPDATE galex SET fuv_lum=%f WHERE gid=%d"%(micJy2SolarLum(fuv,z),gid)
        curs.execute(que)
    conn.commit()
    conn.close()

def fill_Ha_lum():
    conn,curs=setupdb('/home/tom/galex/data.db')
    query=curs.execute("SELECT objID,z,Ha_h,Ha_s FROM sdss")
    for id,z,height,width in query.fetchall():
        que="UPDATE sdss SET Ha_lum=%f WHERE objID=%d"%(sdssflux2Watt(height,width,z)/3.846E26,id) # IN SOLAR LUMINOSITIES
        curs.execute(que)
    conn.commit()
    conn.close()

def fill_agn():
    conn,curs=setupdb('/home/tom/galex/data.db')
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



#
# READ FILES (obsolete)
#
def get_sdss(file='/home/tom/galex/galex2_thomasmarquart.csv'):
    objID,specobjID=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Int64',usecols=[0,1]))
    ra,dec,z,ext_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Float64',usecols=[2,3,4,5,6,7,8,9,10,11]))
    return objID,specobjID,ra,dec,z,ext_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h

def get_galex(file='/home/tom/galex/fluxes_thomasmarquart.csv'):
    sid,gid=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Int64',usecols=[0,1]))
    ra,dec,z,fuv_flux,nuv_flux,fuv_corr=P.transpose(P.load(file,skiprows=1,delimiter=',',dtype='Float64',usecols=[2,3,4,5,6,7]))
    return sid,gid,ra,dec,z,fuv_flux,nuv_flux,fuv_corr
 
