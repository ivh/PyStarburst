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

import Image
import os
from urllib import urlopen
from lxml import etree

DBNAME='lars2.sqlite'

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

def uext2fuv(uext,flux=1.0):
    """ gives the factor to multiply the fuv-flux with, dependign on u-band extiction in magnitudes"""
    return flux*10**(1.66*uext/2.5)

def uext2nuv(uext):
    """ gives the factor to multiply the fuv-flux with, dependign on u-band extiction in magnitudes"""
    return 10**(1.66*uext/2.5)


#
# FILLING CERTAIN DATABASE-TABLES
#

def fill_ext(curs):
    createcolumnifnotexists(curs,'fuv_corr','REAL','galex')
    query=curs.execute("SELECT g.sid,g.gid,s.extinction_u,g.fuv_flux FROM sdss s, galex g WHERE g.sid=s.ObjID")
    sid,gid,ext,fuv = zip(*query.fetchall())
    fc = uext2fuv(N.array(ext),N.array(fuv))
    curs.executemany("UPDATE galex SET fuv_corr=? WHERE gid=?", zip(*(fc,gid)))
    #for sid,gid,ext,fuv in query.fetchall():
    #    que="UPDATE galex SET fuv_corr=%f WHERE gid=%d"%(fuv*uext2fuv(ext),gid)
    #    curs.execute(que)

def fill_fuv_lum(curs):
    createcolumnifnotexists(curs,'fuv_lum','REAL','galex')
    query=curs.execute("SELECT g.sid,g.gid,s.extinction_u,g.fuv_corr,s.z FROM sdss s, galex g WHERE g.sid=s.ObjID")
    sid,gid,ext,fuv,z = zip(*query.fetchall())
    lum=micJy2SolarLum(N.array(fuv),N.array(z))
    curs.executemany("UPDATE galex SET fuv_lum=? WHERE gid=?",zip(*(lum,gid)))
    #for sid,gid,ext,fuv,z in query.fetchall():
    #    que="UPDATE galex SET fuv_lum=%f WHERE gid=%d"%(micJy2SolarLum(fuv,z),gid)
    #    curs.execute(que)

def fill_fuv_int(curs):
    createcolumnifnotexists(curs,'fuv_int','REAL','galex')
    createcolumnifnotexists(curs,'compact','INTEGER','galex')
    query=curs.execute("SELECT g.sid,g.gid,g.fuv_lum,s.z,s.petroR50_u FROM sdss s, galex g WHERE g.sid=s.ObjID")
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
    createcolumnifnotexists(curs,'Ha_lum','REAL')
    query=curs.execute("SELECT sid,z,Ha_h,Ha_s FROM sdss")
    for id,z,height,width in query.fetchall():
        que="UPDATE sdss SET Ha_lum=%f WHERE sid=%d"%(sdssflux2Watt(height,width,z)/3.846E26,id) # IN SOLAR LUMINOSITIES
        curs.execute(que)

def fill_beta(curs):
    createcolumnifnotexists(curs,'beta','REAL','galex')
    query=curs.execute("SELECT gid,fuv_flux,nuv_flux FROM galex")
    c=N.log10(1530.0/2270.0)
    d=(2270.0/1530.0)**2
    for id,fuv,nuv in query.fetchall():
        beta=-1*N.log10(fuv/nuv*d) /c
        if '%f'%beta == 'nan': continue
        #print fuv,nuv,type(beta)
        curs.execute("UPDATE galex SET beta=%f WHERE gid=%d"%(beta,id))

def fill_agn(curs):
    createcolumnifnotexists(curs,'agn','INTEGER')
    ids=gettable(curs,cols='objID',where='(Ha_h >0) AND (Hb_h>0) AND (OIII_h>0) AND (NII_h>0)',table='sdss')
    x,y,sig=gettable(curs,cols='NII_h/Ha_h,OIII_h/Hb_h,Ha_s',where='(Ha_h >0) AND (Hb_h>0) AND (OIII_h>0) AND (NII_h>0)',table='sdss')
    x=N.log10(N.array(x))
    y=N.log10(N.array(y))
    agn=N.where( (y>mylee) | (sig > 5),1,0)
    for i,id in enumerate(ids):
        curs.execute("UPDATE sdss SET agn=%s WHERE objID=%s"%(agn[i],id))

def fill_from_csv(curs, filename, tablename, nint=2):
    """
    create columns according to file header and fill them
    """
    f = open(filename)
    cols = f.readline().split(',')
    types = ['INTEGER']*nint + ['REAL']*(len(cols)-nint)
    coldef = ','.join(['%s %s'%(col,types[i]) for i,col in enumerate(cols)])

    curs.execute('CREATE TABLE IF NOT EXISTS %s (%s);'%(tablename,coldef))

    for line in f:
        curs.execute('INSERT INTO %s VALUES (%s)'%(tablename,','.join(['?']*len(cols))), line.split(','))

def fill_sdss_name(curs,table='sel'):
    createcolumnifnotexists(curs,'name','TEXT',table=table)
    urlbase='http://cas.sdss.org/dr7/en/tools/explore/summary.asp?spec=&id='
    for sid,z in curs.execute('select objID,z from %s'%table).fetchall():
        url=urlbase+hex(int(sid))
        u = urlopen(url)
        h=etree.parse(u,parser=etree.HTMLParser())
        for i in h.iter():
            if i.tag=='h2':
                for j in i.getchildren():
                    if j.text: name=j.text.split()[-1]
        curs.execute('update %s SET name="%s" where objID=%s'%(table,name,sid))


def makedb(dbname=DBNAME,sdss='sdss.csv',sints=2,galex='galex.csv',gints=3):
    conn,curs=setupdb(dbname)
    fill_from_csv(curs,sdss,'sdss',nint=sints)
    fill_from_csv(curs,galex,'galex',nint=gints)

    curs.execute('delete from sdss where extinction_u IS NULL;')
    curs.execute('delete from galex where specid not in (select specobjid from sdss);')
    conn.commit()

    # do the work
    curs.execute('BEGIN TRANSACTION;')
    fill_ext(curs)
    fill_fuv_lum(curs)
    curs.execute('delete from galex where fuv_lum ISNULL;')
    curs.execute('delete from galex where fuv_s2n = -999.0;')
    curs.execute('delete from galex where fuv_lum<1E10 and specid in (select specObjID from sdss where z > 0.05)')
    curs.execute('delete from sdss where specObjID not in (select specid from galex )')
    conn.commit()
    fill_fuv_int(curs)
    fill_agn(curs)
    fill_beta(curs)

    curs.execute('vacuum;')
    conn.commit()
    return conn

def selection2html(curs,outfile='sel2013.html',where='',table='sel'):
    urlbase='http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id='
    wanted='name,objID, gid, ra, dec, z, fuv_lum, fuv_s2n,fuv_int, extinction_u, beta, agn, Ha_w, Ha_s, Ha_h,Hb_h, Hd_w, OIII_h,NII_h,mag_g,mag_r'
    wanto=wanted.replace('extinction','A')
    wants=wanto.split(',')
    f=open(outfile,'w')
    f.write('<table border=1><tr><th>#</th>\n')
    for w in wants: f.write('<th>%s</th>'%w)
    f.write('<th>M_g</th><th>M_r</th>')
    f.write('</tr>\n')
    if where: where='AND %s'%where
    curs.execute('SELECT DISTINCT %s FROM %s %s ORDER BY fuv_lum desc'%(wanted,table,where))
    data=curs.fetchall()
    prevID,counter=0,1
    for name,sid,gid,ra,dec,z,fuv_lum,fuv_s2n,fuv_int,A_u,beta,agn,Ha_w,Ha_s,Ha_h,Hb_h,Hd_w,OIII_h,NII_h,mag_g,mag_r in data:
        f.write('<tr><td>%s</td><td>%s</td>'%(counter,name))
        if sid!=prevID: f.write('<td><a href="%s">%s</a></td>'%(urlbase+str(sid),str(sid)))
        else: f.write('<td></td>')
        f.write('<td>%s</td><td>%.5f</td><td>%.5f</td><td>%.4f</td>'%(gid,ra,dec,z))
        f.write('<td>%.3f</td><td>%.1f</td><td>%.3f</td>'%(N.log10(fuv_lum),fuv_s2n,N.log10(fuv_int)))
        f.write('<td>%.2f</td><td>%.2f</td><td>%s</td>'%(A_u,beta or N.nan,agn))
        f.write('<td>%.1f</td><td>%.1f</td><td>%.1f</td>'%(Ha_w,Ha_s,Ha_h))
        f.write('<td>%.1f</td><td>%.1f</td><td>%.1f</td><td>%.1f</td>'%(Hb_h,Hd_w,OIII_h,NII_h))
        f.write('<td>%.1f</td><td>%.1f</td>'%(mag_g,mag_r))
        f.write('<td>%.1f</td>'%(absmag(mag_r,z)))
        f.write('<td>%.1f</td>'%(absmag(mag_g,z)))

        f.write('</tr>\n')
        prevID=sid
        counter+=1
    f.close()

def selection2ascii(curs,outfile='sel2013.dat',where=''):
    wanted='s.ObjID, g.gid, s.ra, s.dec, s.z, g.fuv_lum, g.fuv_s2n,g.fuv_int, s.extinction_u, g.beta, s.agn, s.Ha_w, s.Ha_s, s.Ha_h,s.Hb_h, Hd_w, s.OIII_h,s.NII_h'
    wanto=wanted.replace('s.','').replace('g.','').replace('extinction','A')
    f=open(outfile,'w')
    f.write('# '+ wanto)
    if where !='': where='AND %s'%where
    curs.execute('SELECT DISTINCT %s FROM sdss s, galex g WHERE g.sid=s.ObjID %s ORDER BY s.z'%(wanted,where))
    data=curs.fetchall()
    for sid,gid,ra,dec,z,fuv_lum,fuv_s2n,fuv_int,A_u,beta,agn,Ha_w,Ha_s,Ha_h,Hb_h,Hd_w,OIII_h,NII_h in data:
        f.write('%s, %s, '%(sid,sid))
        f.write('%s, %.3f, %.3f, %.3f, '%(gid,ra,dec,z))
        f.write('%.3f, %.1f, %.3f, '%(N.log10(fuv_lum),fuv_s2n,N.log10(fuv_int)))
        f.write('%.2f, %.2f, %s, '%(A_u,beta or N.nan,agn))
        f.write('%.1f, %.1f, %.1f, '%(Ha_w,Ha_s,Ha_h))
        f.write('%.1f, %.1f, %.1f, %.1f\n'%(Hb_h,Hd_w,OIII_h,NII_h))
    f.close()


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

    # plot the dividing line
    xl=P.arange(-1.7,-0.15,0.01)
    yl=sdss.mylee(xl)
    P.plot(xl,yl,'r-',linewidth=2)
    P.xlabel('log [NII]/Ha')
    P.ylabel('log [OIII]/Hb')

def xyZFUV(curs,where):
    return gettable(curs,'s.z,g.fuv_lum',table='sdss s, galex g',where='g.sid=s.objID AND %s'%where)

def plot_z_fuv(curs):
    z,fuv=xyZFUV(curs,'s.agn=0')
    P.plot(z,N.log10(fuv),'k,')
    #z,fuv=xyZFUV(curs,'g.compact=1')
    #P.semilogy(z,fuv,'y,')
    P.xlabel('z')
    P.ylabel(r'L$_{FUV}$')
    #P.legend(loc='lower right')

def plot_z_fuv_prop(curs):
    z,fuv=xyZFUV(curs,'s.agn=0')
    P.plot(z,N.log10(fuv),'k,',alpha=0.1)
    P.xlabel('$z$')
    P.ylabel(r'$\log(L_{FUV})$')

def xyHaLum(curs,where):
    return gettable(curs,'s.Ha_w,g.fuv_lum',table='sdss s, galex g',where='g.sid=s.objID AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_Ha_lum(curs):
    Ha,fuv=xyHaLum(curs,'g.compact=0')
    P.semilogy(Ha,fuv,'b,')
    Ha,fuv=xyHaLum(curs,'g.compact=1')
    P.loglog(Ha,fuv,'y,')
    Ha,fuv=xyHaLum(curs,'g.compact=2')
    P.loglog(Ha,fuv,'r,')
    P.xlabel(r'$\mathrm{EW(H\alpha)}$')
    P.ylabel(r'L$_{FUV}$')
    P.legend(loc='lower right')


def xyHaInt(curs,where):
    return gettable(curs,'s.Ha_w,g.fuv_int',table='sdss s, galex g',where='g.sid=s.objID AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_Ha_int(curs):
    Ha,fuv=xyHaInt(curs,'g.compact=0')
    P.loglog(Ha,fuv,'b,')
    Ha,fuv=xyHaInt(curs,'g.compact=1')
    P.loglog(Ha,fuv,'y,')
    Ha,fuv=xyHaInt(curs,'g.compact=2')
    P.loglog(Ha,fuv,'r,')
    P.xlabel(r'$\mathrm{EW(H\alpha)}$')
    P.ylabel(r'I$_{FUV}$')
    P.legend(loc='lower right')

def xyHaBeta(curs,where):
    return gettable(curs,'s.Ha_w,g.beta',table='sdss s, galex g',where='g.sid=s.objID AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_Ha_beta(curs):
    Ha,beta=xyHaBeta(curs,'g.compact=0')
    P.semilogx(Ha,beta,'b,')
    Ha,beta=xyHaBeta(curs,'g.compact=1')
    P.semilogx(Ha,beta,'y,')
    Ha,beta=xyHaBeta(curs,'g.compact=2')
    P.semilogx(Ha,beta,'r,')
    P.xlabel(r'$\mathrm{EW(H\alpha)}$')
    P.ylabel(r'$\beta$')

def xylumBeta(curs,where):
    return gettable(curs,'g.fuv_lum,g.beta',table='sdss s, galex g',where='g.sid=s.objID AND g.fuv_corr NOT NULL AND s.z NOT NULL AND s.agn=0 AND %s'%where)
def plot_lum_beta(curs):
    lum,beta=xylumBeta(curs,'g.compact=0')
    P.semilogx(lum,beta,'b,')
    lum,beta=xylumBeta(curs,'g.compact=1')
    P.semilogx(lum,beta,'y,')
    lum,beta=xylumBeta(curs,'g.compact=2')
    P.semilogx(lum,beta,'r,')
    P.xlabel(r'L$_{FUV}$')
    P.ylabel(r'$\beta$')

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

    P.subplots_adjust(0.07,.11,.98,.98,.27,.24)


def getselimages(curs,table='sel'):
    url='http://casjobs.sdss.org/ImgCutoutDR7/getjpeg.aspx?ra=%s&dec=%s&scale=0.19806&width=256&height=256'

    curs.execute('SELECT objID,ra,dec from %s'%table)
    data=curs.fetchall()
    for sid,ra,dec in data:
        fname='obj_%s.jpg'%str(sid)
        if os.path.exists(fname):
            continue
        im=urlopen(url%(ra,dec))
        f=open(fname,'w')
        f.write(im.read())
        f.close()
        im.close()

def fooplot(curs,ax1=None,ax2=None):
    #P.figure(figsize=(10,4.96))
    if not ax1: ax1=P.axes([0.38,0.01,0.30,0.32])
    if not ax2: ax2=P.axes([0.68,0.01,0.31,0.32],sharey=ax1)
    P.rc('xtick.major', pad=-12)
    P.rc('xtick', labelsize=10)

    curs.execute('SELECT DISTINCT sid,z,Ha_w from sel ORDER BY fuv_lum desc')
    data=curs.fetchall()
    for sid,z,ha in data:
        fuv,beta=curs.execute('select max(fuv_lum),beta from sel where objID=%s'%sid).fetchone()

        if ha < 70: color='g';marker='o'
        elif ha > 70 and ha < 140: color='y';marker='s'
        elif ha > 140: color='r';marker='D'
        ax1.plot(-1*(beta or N.nan),N.log10(fuv),color=color,marker=marker)
        ax2.plot(N.log10(ha),N.log10(fuv),color=color,marker=marker)

    P.setp(ax1.get_yticklabels(), visible=False)
    P.setp(ax2.get_yticklabels(), fontsize=10)
    ax1.xaxis.labelpad=-19
    ax2.xaxis.labelpad=-19
    ax2.set_xlabel(r'$\log\,W(H_\alpha)$')
    ax1.set_ylabel(r'$\log(L_{FUV})$')
    ax1.set_xlabel(r'$\beta$')
    ax1.axis([-0.4,-1.8,8.75,10.45])
    ax2.axis([1.55,2.45,8.75,10.45])
    ax1.set_yticks([8.8, 9.2, 9.6,10,10.4])
    ax1.set_xticks([-0.5,-1.0,-1.5])
    ax2.set_xticks([1.8,2.0,2.2,2.4])

def plotselimages(curs):
    P.figure(figsize=(8.66,7.5))
    P.subplots_adjust(0.01,0.01,0.99,0.99,0.02,0.02)
    curs.execute('SELECT DISTINCT sid,z,Ha_w from sel ORDER BY fuv_lum desc')
    data=curs.fetchall()
    i=1
    for sid,z,ha in data:
        ax=P.subplot(6,7,i,frameon=False)
        ax.set_axis_off()
        P.imshow(Image.open('obj_%s.jpg'%str(sid)),origin='lower')
        fuv=curs.execute('select max(fuv_lum) from sel where objID=%s'%sid).fetchone()
        P.text(6,3,'$\mathbf{%d}$'%i,fontsize=18,color='w')
        P.text(100,7,'$z:\\, %.3f$'%z,fontsize=11,color='w')
        P.text(5,215,'$W(H\\alpha):\\, %d$'%ha,fontsize=11,color='w')
        P.text(5,185,'$\log(L_{FUV}):\\, %.1f $'%N.log10(fuv),fontsize=11,color='w')
        i+=1

    ax1=P.axes([0.18,0.01,0.40,0.3])
    ax2=P.axes([0.58,0.01,0.40,0.3],sharey=ax1)

    fooplot(curs,ax1,ax2)

