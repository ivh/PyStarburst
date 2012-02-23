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

DBNAME='data.db'

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

def uext2nuv(uext):
    """ gives the factor to multiply the fuv-flux with, dependign on u-band extiction in magnitudes"""
    return 10**(1.66*uext/2.5)


#
# PLOTTING
#
def xyAGN(curs,where):
    Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h=gettable(curs,'Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h',table='sdss',where=where)
    x=P.log10(NII_h/Ha_h)
    y=P.log10(OIII_h/Hb_h)
    return x,y
def plotAGN(curs):
    x,y=xyAGN(curs,where='agn=1')
    P.plot(x,y,'r,')
    x,y=xyAGN(curs,where='agn=0')
    P.plot(x,y,'b,')

    x,y=xyAGN(curs,where='sid in (SELECT sid from mccand)')
    P.plot(x,y,'ob')
    x,y=xyAGN(curs,where='sid in (SELECT sid from green)')
    P.plot(x,y,'*g')
    x,y=xyAGN(curs,where='sid in (SELECT sid from heck3)')
    P.plot(x,y,'sr')
    x,y=xyAGN(curs,where='sid in (SELECT sid from heck2)')
    P.plot(x,y,'Dr')
    x,y=xyAGN(curs,where='sid in (SELECT sid from heck1)')
    P.plot(x,y,'*r')

    # plot the dividing line
    xl=P.arange(-1.7,-0.15,0.01)
    yl=sdss.mylee(xl)
    P.plot(xl,yl,'r-',linewidth=2)
    P.xlabel('log [NII]/Ha')
    P.ylabel('log [OIII]/Hb')

def xyZFUV(curs,where):
    return gettable(curs,'s.z,g.fuv_lum',table='sdss s, galex g',where='g.sid=s.sid AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_z_fuv(curs):
    z,fuv=xyZFUV(curs,'g.compact=0')
    P.semilogy(z,fuv,'b,')
    z,fuv=xyZFUV(curs,'g.compact=1')
    P.semilogy(z,fuv,'y,')
    z,fuv=xyZFUV(curs,'g.compact=2')
    P.semilogy(z,fuv,'r,')
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from mccand)')
    P.semilogy(z,fuv,'ob')
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from green)')
    P.semilogy(z,fuv,'*g')
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from heck3)')
    P.semilogy(z,fuv,'sr')
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from heck2)')
    P.semilogy(z,fuv,'Dr')
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from heck1)')
    P.semilogy(z,fuv,'*r')
    P.xlabel('z')
    P.ylabel(r'L$_{FUV}$')
    #P.legend(loc='lower right')

def plot_z_fuv_prop(curs):
    z,fuv=xyZFUV(curs,'s.agn=0')
    P.plot(z,N.log10(fuv),'k,',alpha=0.1)
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from sel)')
    P.plot(z,N.log10(fuv),'ob')
    z,fuv=xyZFUV(curs,'s.sid IN (SELECT sid from sel) AND (s.sid IN (SELECT sid from heck3) OR s.sid IN (SELECT sid from heck2) OR s.sid IN (SELECT sid from heck1) OR s.sid IN (SELECT sid from green) OR s.sid IN (SELECT sid from mccand))')
    P.plot(z,N.log10(fuv),'or')
    P.xlabel('$z$')
    P.ylabel(r'$\log(L_{FUV})$')

def xyHaLum(curs,where):
    return gettable(curs,'s.Ha_w,g.fuv_lum',table='sdss s, galex g',where='g.sid=s.sid AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_Ha_lum(curs):
    Ha,fuv=xyHaLum(curs,'g.compact=0')
    P.semilogy(Ha,fuv,'b,')
    Ha,fuv=xyHaLum(curs,'g.compact=1')
    P.loglog(Ha,fuv,'y,')
    Ha,fuv=xyHaLum(curs,'g.compact=2')
    P.loglog(Ha,fuv,'r,')
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from mccand)')
    P.loglog(Ha,fuv,'bo')
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from green)')
    P.loglog(Ha,fuv,'*g')
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from heck3)')
    P.loglog(Ha,fuv,'sr')
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from heck2)')
    P.loglog(Ha,fuv,'Dr')
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from heck1)')
    P.loglog(Ha,fuv,'*r')
    P.xlabel('Ha_W')
    P.ylabel(r'L$_{FUV}$')
    P.legend(loc='lower right')

def plot_LyHa_prop(curs):
    Ha,fuv=xyHaLum(curs,'s.agn=0')
    P.plot(N.log10(Ha),N.log10(fuv),',k',alpha=0.1)
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from sel)')
    P.plot(N.log10(Ha),N.log10(fuv),'bo')
    Ha,fuv=xyHaLum(curs,'s.sid IN (SELECT sid from sel) AND (s.sid IN (SELECT sid from heck3) OR s.sid IN (SELECT sid from heck2) OR s.sid IN (SELECT sid from heck1) OR s.sid IN (SELECT sid from green) OR s.sid IN (SELECT sid from mccand))')
    P.plot(N.log10(Ha),N.log10(fuv),'ro')
    P.xlabel(r'$\log\,W(H_\alpha)$')
#P.ylabel(r'$L_{FUV}$')


def xyHaInt(curs,where):
    return gettable(curs,'s.Ha_w,g.fuv_int',table='sdss s, galex g',where='g.sid=s.sid AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_Ha_int(curs):
    Ha,fuv=xyHaInt(curs,'g.compact=0')
    P.loglog(Ha,fuv,'b,')
    Ha,fuv=xyHaInt(curs,'g.compact=1')
    P.loglog(Ha,fuv,'y,')
    Ha,fuv=xyHaInt(curs,'g.compact=2')
    P.loglog(Ha,fuv,'r,')
    Ha,fuv=xyHaInt(curs,'s.sid IN (SELECT sid from mccand)')
    P.loglog(Ha,fuv,'bo')
    Ha,fuv=xyHaInt(curs,'s.sid IN (SELECT sid from green)')
    P.loglog(Ha,fuv,'*g')
    Ha,fuv=xyHaInt(curs,'s.sid IN (SELECT sid from heck3)')
    P.loglog(Ha,fuv,'sr')
    Ha,fuv=xyHaInt(curs,'s.sid IN (SELECT sid from heck2)')
    P.loglog(Ha,fuv,'Dr')
    Ha,fuv=xyHaInt(curs,'s.sid IN (SELECT sid from heck1)')
    P.loglog(Ha,fuv,'*r')
    P.xlabel(r'W(H$_\alpha$)')
    P.ylabel(r'I$_{FUV}$')
    P.legend(loc='lower right')

def xyHaBeta(curs,where):
    return gettable(curs,'s.Ha_w,g.beta',table='sdss s, galex g',where='g.sid=s.sid AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_Ha_beta(curs):
    Ha,beta=xyHaBeta(curs,'g.compact=0')
    P.semilogx(Ha,beta,'b,')
    Ha,beta=xyHaBeta(curs,'g.compact=1')
    P.semilogx(Ha,beta,'y,')
    Ha,beta=xyHaBeta(curs,'g.compact=2')
    P.semilogx(Ha,beta,'r,')
    Ha,beta=xyHaBeta(curs,'s.sid IN (SELECT sid from mccand)')
    P.semilogx(Ha,beta,'bo')
    Ha,beta=xyHaBeta(curs,'s.sid IN (SELECT sid from green)')
    P.semilogx(Ha,beta,'*g')
    Ha,beta=xyHaBeta(curs,'s.sid IN (SELECT sid from heck3)')
    P.semilogx(Ha,beta,'sr')
    Ha,beta=xyHaBeta(curs,'s.sid IN (SELECT sid from heck2)')
    P.semilogx(Ha,beta,'Dr')
    Ha,beta=xyHaBeta(curs,'s.sid IN (SELECT sid from heck1)')
    P.semilogx(Ha,beta,'*r')
    P.xlabel('Ha_W')
    P.ylabel('beta')

def xylumBeta(curs,where):
    return gettable(curs,'g.fuv_lum,g.beta',table='sdss s, galex g',where='g.sid=s.sid AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_lum_beta(curs):
    lum,beta=xylumBeta(curs,'g.compact=0')
    P.semilogx(lum,beta,'b,')
    lum,beta=xylumBeta(curs,'g.compact=1')
    P.semilogx(lum,beta,'y,')
    lum,beta=xylumBeta(curs,'g.compact=2')
    P.semilogx(lum,beta,'r,')
    lum,beta=xylumBeta(curs,'s.sid IN (SELECT sid from mccand)')
    P.semilogx(lum,beta,'bo')
    lum,beta=xylumBeta(curs,'s.sid IN (SELECT sid from green)')
    P.semilogx(lum,beta,'*g')
    lum,beta=xylumBeta(curs,'s.sid IN (SELECT sid from heck3)')
    P.semilogx(lum,beta,'sr')
    lum,beta=xylumBeta(curs,'s.sid IN (SELECT sid from heck2)')
    P.semilogx(lum,beta,'Dr')
    lum,beta=xylumBeta(curs,'s.sid IN (SELECT sid from heck1)')
    P.semilogx(lum,beta,'*r')
    P.xlabel(r'L$_{FUV}$')
    P.ylabel('beta')

def plotprop(curs):
    P.figure(figsize=(10,4))
    P.subplots_adjust(0.09,0.15,0.97,0.95,0.0,0.0)
    ax1=P.subplot(121)
    plot_z_fuv_prop(curs)
    ax2=P.subplot(122,sharey=ax1)
    plot_LyHa_prop(curs)
    P.setp(ax2.get_yticklabels(), visible=False)
    ax1.axis([0.025,0.185,8.95,10.78])
    ax2.axis([1.9,3.01,8.95,10.78])

def plotall(curs):
    P.clf()
    P.subplot(231)
    plotAGN(curs)
    P.subplot(232)
    plot_z_fuv(curs)
    P.subplot(233)
    plot_Ha_lum(curs)
    P.subplot(234)
    plot_Ha_int(curs)
    P.subplot(235)
    plot_Ha_beta(curs)
    P.subplot(236)
    plot_lum_beta(curs)

def detcolor(curs,sid):
    curs.execute('select sid from green where sid="%d"'%sid)
    if curs.fetchone(): return '#48db00'
    curs.execute('select sid from mccand where sid="%d"'%sid)
    if curs.fetchone(): return 'r'
    curs.execute('select sid from heck where sid="%d"'%sid)
    if curs.fetchone(): return '#429fff'

    return 'y'


def plotselimages(curs):
    import Image
    P.rc('xtick.major', pad=-12)
    P.rc('xtick', labelsize=8)
    P.figure(figsize=(10,4.96))
    P.subplots_adjust(0.01,0.01,0.99,0.99,0.02,0.02)

    ax1=P.axes([0.38,0.01,0.30,0.32])
    #plot_z_fuv_prop(curs)
    ax2=P.axes([0.68,0.01,0.31,0.32],sharey=ax1)
    #plot_LyHa_prop(curs)
    P.setp(ax1.get_yticklabels(), visible=False)
    P.setp(ax2.get_yticklabels(), fontsize=8)
    ax1.set_yticks([9.5,10,10.5])
    ax1.set_xticks([0.04,0.08,0.12,0.16])
    ax2.set_xticks([2,2.4,2.8])
    ax1.xaxis.labelpad=-19
    ax2.xaxis.labelpad=-19
    ax2.set_xlabel(r'$\log\,W(H_\alpha)$')
    ax1.set_ylabel(r'$\log(L_{FUV})$')
    ax1.set_xlabel(r'$z$')

    curs.execute('SELECT DISTINCT sid,z,Ha_w from sdss WHERE (sid IN (SELECT sid FROM sel)) ORDER BY z')
    data=curs.fetchall()
    i=1
    for sid,z,ha in data:
        ax=P.subplot(3,6,i,frameon=False)
        ax.set_axis_off()
        P.imshow(Image.open('%s.jpg'%str(sid)),origin='lower')
        color=detcolor(curs,sid)
        P.text(6,3,'$\mathbf{%d}$'%i,fontsize=18,color=color)
        P.text(165,7,'$z:\\, %.3f$'%z,fontsize=11,color='w')
        P.text(5,215,'$W(H_\\alpha):\\, %d \\AA$'%ha,fontsize=11,color='w')
        curs.execute('SELECT fuv_lum from galex WHERE sid=%s'%str(sid))
        fuv=N.max(curs.fetchall())
        P.text(5,185,'$\log(L_{FUV}):\\, %.1f $'%N.log10(fuv),fontsize=11,color='w')
        ax1.plot(z,N.log10(fuv),color=color,marker='o')
        ax2.plot(N.log10(ha),N.log10(fuv),color=color,marker='o')
        i+=1

    ax1.axis([0.02,0.185,9.0,10.78])
    ax2.axis([1.9,2.9,9.0,10.78])

def getselimages(curs):
    url='http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx?ra=%s&dec=%s&scale=0.19806&width=256&height=256'

    from urllib import urlopen as get
    curs.execute('SELECT sid,ra,dec from sdss WHERE sid IN (SELECT sid FROM sel)')
    data=curs.fetchall()
    for sid,ra,dec in data:
        im=get(url%(ra,dec))
        f=open('%s.jpg'%str(sid),'w')
        f.write(im.read())
        f.close()
        im.close()


#
# FILLING CERTAIN DATABASE-TABLES
#

def fill_ext(curs):
    query=curs.execute("SELECT g.sid,g.gid,s.ext_u,g.fuv FROM sdss s, galex g WHERE g.sid=s.sid")
    for sid,gid,ext,fuv in query.fetchall():
        que="UPDATE galex SET fuv_corr=%f WHERE gid=%d"%(fuv*uext2fuv(ext),gid)
        curs.execute(que)

def fill_fuv_lum(curs):
    query=curs.execute("SELECT g.sid,g.gid,s.ext_u,g.fuv_corr,s.z FROM sdss s, galex g WHERE g.sid=s.sid")
    for sid,gid,ext,fuv,z in query.fetchall():
        que="UPDATE galex SET fuv_lum=%f WHERE gid=%d"%(micJy2SolarLum(fuv,z),gid)
        curs.execute(que)

def fill_fuv_int(curs):
    query=curs.execute("SELECT g.sid,g.gid,g.fuv_lum,s.z,s.petroR50_u FROM sdss s, galex g WHERE g.sid=s.sid")
    for sid,gid,fuv_lum,z,r in query.fetchall():
        fuv_int=fuv_lum/ (N.pi * (r**2))
        curs.execute("UPDATE galex SET fuv_int=%f WHERE gid=%d"%(fuv_int,gid))
        if fuv_int>10**9:
            curs.execute("UPDATE galex SET compact=2 WHERE gid=%d"%gid)
        elif fuv_int>10**8:
            curs.execute("UPDATE galex SET compact=1 WHERE gid=%d"%gid)
        else:
            curs.execute("UPDATE galex SET compact=0 WHERE gid=%d"%gid)

def fill_Ha_lum(curs):
    query=curs.execute("SELECT sid,z,Ha_h,Ha_s FROM sdss")
    for id,z,height,width in query.fetchall():
        que="UPDATE sdss SET Ha_lum=%f WHERE sid=%d"%(sdssflux2Watt(height,width,z)/3.846E26,id) # IN SOLAR LUMINOSITIES
        curs.execute(que)

def fill_beta(curs):
    query=curs.execute("SELECT gid,fuv,nuv FROM galex")
    c=N.log10(1530.0/2270.0)
    d=(2270.0/1530.0)**2
    for id,fuv,nuv in query.fetchall():
        beta=-1*N.log10(fuv/nuv*d) /c
        #print fuv,nuv,beta
        curs.execute("UPDATE galex SET beta=%f WHERE gid=%d"%(beta,id))

def fill_agn(curs):
    ids=gettable(curs,cols='sid',where='(Ha_h >0) AND (Hb_h>0) AND (OIII_h>0) AND (NII_h>0)',table='sdss')
    x,y,sig=gettable(curs,cols='NII_h/Ha_h,OIII_h/Hb_h,Ha_s',where='(Ha_h >0) AND (Hb_h>0) AND (OIII_h>0) AND (NII_h>0)',table='sdss')
    x=N.log10(N.array(x))
    y=N.log10(N.array(y))
    agn=N.where( (y>mylee) | (sig > 5),1,0)
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET agn=%s WHERE sid=%s"%(agn[i],id))


def fill_galex(curs,galex):
    """
    galex columns:sid,gid,distance,xnum,fuv,nuv

    """
    sid,gid,xnum=N.loadtxt(galex,skiprows=1,unpack=True,dtype='S',delimiter=',',usecols=(0,1,3))
    distance,fuv,nuv=N.loadtxt(galex,skiprows=1,unpack=True,dtype='Float64',delimiter=',',usecols=(2,4,5))
    curs.execute('CREATE TABLE IF NOT EXISTS galex (gid INTEGER, sid INTEGER, fuv REAL, nuv REAL, xnum INTEGER, distance REAL, fuv_corr REAL, fuv_lum REAL, fuv_int REAL, compact INTEGER, beta REAL);')
    for i,id in enumerate(gid):
        d=(id,sid[i],fuv[i],nuv[i],xnum[i],distance[i])
        if (fuv[i]>0) and (nuv[i]>0):
            curs.execute('INSERT INTO galex VALUES (%s,%s,%f,%f,%s,%f,NULL, NULL, NULL, NULL, NULL);'%d)

def fill_sdss(curs,sdss):
    """
    sdss columns: objID,specobjID,ra,dec,l,b,z,mag_u,mag_g,mag_r,mag_i,mag_z,petroRad_u,petroR50_u,extinction_u,isoA_u,isoB_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h
    """
    objID,specobjID=N.loadtxt(sdss,skiprows=1,unpack=True,dtype='S',delimiter=',',usecols=(0,1))
    ra,dec,l,b,z,mag_u,mag_g,mag_r,mag_i,mag_z,petroRad_u,petroR50_u,extinction_u,isoA_u,isoB_u,Ha_w,Ha_h,Ha_s,Hb_h,OIII_h,NII_h=N.loadtxt(sdss,skiprows=1,unpack=True,dtype='Float64',delimiter=',',usecols=tuple(N.arange(21)+2))
    curs.execute('CREATE TABLE IF NOT EXISTS sdss (sid INTEGER, specid INTEGER, ra REAL, dec REAL, l REAL, b REAL, z REAL, mag_u REAL, mag_g REAL, mag_r REAL, mag_i REAL, mag_z REAL, petroRad_u REAL, petroR50_u REAL, ext_u REAL, isoA_u REAL, isoB_u REAL, Ha_w REAL, Ha_h REAL, Ha_s REAL, Hb_h REAL, OIII_h REAL, NII_h REAL, agn INTEGER DEFAULT 0, Ha_lum REAL);')
    for i,id in enumerate(objID):
        d=(id,specobjID[i],ra[i],dec[i],l[i],b[i],z[i],mag_u[i],mag_g[i],mag_r[i],mag_i[i],mag_z[i],petroRad_u[i],petroR50_u[i],extinction_u[i],isoA_u[i],isoB_u[i],Ha_w[i],Ha_h[i],Ha_s[i],Hb_h[i],OIII_h[i],NII_h[i])
        curs.execute('INSERT INTO sdss VALUES (%s,%s,%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, NULL, NULL);'%d)

def makedb(dbname=DBNAME,sdss='sdss.csv',galex='galex.csv'):
    conn,curs=setupdb(dbname)
    fill_galex(curs,galex)
    conn.commit()
    fill_sdss(curs,sdss)
    conn.commit()

    #fill_galex(curs,'galex_foreign.csv')
    #fill_sdss(curs,'sdss_foreign.csv')
    #conn.commit()

    # get the "external" runs in..
    #for run in ['green','mccand','heck1','heck2','heck3','sel']:
    #    curs.execute('CREATE TABLE %s (sid INTEGER);'%run)
    #    sid=N.loadtxt('%s_targets'%run,skiprows=1,unpack=True,dtype='S',usecols=(0,))
    #    for id in sid:
    #        curs.execute('INSERT INTO %s VALUES (%s)'%(run,id))
    #curs.execute('create view heck as select * from heck1 union select * from heck2 union select * from heck3')
    #conn.commit()

    return conn
    # do the work
    fill_agn(curs)
    fill_beta(curs)
    fill_ext(curs)
    fill_fuv_lum(curs)
    fill_fuv_int(curs)

    conn.commit()
    return conn

def dump_selection(curs,outfile='selection.html',where=''):
    urlbase='http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id='
    wanted='s.sid, s.z, g.fuv_lum, g.fuv_int, g.beta, s.agn, s.Ha_w, s.Ha_s, s.Ha_h,s.Hb_h,s.OIII_h,s.NII_h'
    wanto=wanted.replace('s.','').replace('g.','')
    wants=wanto.split(',')
    f=open(outfile,'w')
    f.write('<table border=1>\n')
    for w in wants[:-4]: f.write('<td>%s</td>'%w)
    f.write('<td>log(NII/Ha)</td><td>log(OIII/Hb)</td></tr>')
    if where !='': where='AND %s'%where
    curs.execute('SELECT DISTINCT %s FROM sdss s, galex g WHERE g.sid=s.sid %s ORDER BY s.z'%(wanted,where))
    data=curs.fetchall()
    for sid,z,fuv_lum,fuv_int,beta,agn,Ha_w,Ha_s,Ha_h,Hb_h,OIII_h,NII_h in data:
        f.write('<tr>')
        f.write('<td><a href="%s">%s</a></td>'%(urlbase+str(sid),str(sid)))
        f.write('<td>%.3f</td><td>%.3f</td><td>%.3f</td>'%(z,N.log10(fuv_lum),N.log10(fuv_int)))
        f.write('<td>%.2f</td><td>%s</td>'%(beta,str(agn)))
        f.write('<td>%.1f</td><td>%.1f</td>'%(Ha_w,Ha_s))
        x,y=P.log10(NII_h/Ha_h),P.log10(OIII_h/Hb_h)
        f.write('<td>%.2f</td><td>%.2f</td>'%(x,y))

        f.write('</tr>')
    f.close()


def dump_foreign(curs,outfile='foreign.html'):
    urlbase='http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id='
    wanted='s.sid, s.z, g.fuv_lum, g.fuv_int, g.beta, s.agn, s.Ha_w, s.Ha_s, s.Ha_h,s.Hb_h,s.OIII_h,s.NII_h'
    wanto=wanted.replace('s.','').replace('g.','')
    wants=wanto.split(',')
    f=open(outfile,'w')
    f.write('<table border=1><th>\n')
    for run in ['green','mccand','heck1','heck2','heck3']:
        f.write('<tr><td colspan=6><strong>%s</strong></td></tr><tr>'%run)
        for w in wants[:-4]: f.write('<td>%s</td>'%w)
        f.write('<td>log(NII/Ha)</td><td>log(OIII/Hb)</td></tr>')
        curs.execute('SELECT DISTINCT %s FROM sdss s, galex g WHERE g.sid=s.sid AND (s.sid in (SELECT sid from %s)) ORDER BY s.sid'%(wanted,run))
        data=curs.fetchall()
        for sid,z,fuv_lum,fuv_int,beta,agn,Ha_w,Ha_s,Ha_h,Hb_h,OIII_h,NII_h in data:
            f.write('<tr>')
            f.write('<td><a href="%s">%s</a></td>'%(urlbase+str(sid),str(sid)))
            f.write('<td>%.3f</td><td>%.3f</td><td>%.3f</td>'%(z,N.log10(fuv_lum),N.log10(fuv_int)))
            f.write('<td>%.2f</td><td>%s</td>'%(beta,str(agn)))
            f.write('<td>%.1f</td><td>%.1f</td>'%(Ha_w,Ha_s))
            x,y=P.log10(NII_h/Ha_h),P.log10(OIII_h/Hb_h)
            f.write('<td>%.2f</td><td>%.2f</td>'%(x,y))

            f.write('</tr>')
        f.write('<tr><td colspan=6>&nbsp;</td></tr>')
    f.close()

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

