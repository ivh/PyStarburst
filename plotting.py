#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import numpy as N
import pylab as P

from db import *
import sdss
from PyGalKin.tool import smooth_gauss


def plotspecwhere(cursor,where,table='sdss'):
    ids=cursor.execute('SELECT objID from %s WHERE %s'%(table,where))
    for id in ids:
       plotspecbyid(id) 

       
def plotspecbyid(id):
    fits=specfromid(id)
    head,spec,noise=sdss.splitfits(fits)
    wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
    P.plot(wave,spec,linestyle='steps')

def plotfromdb(cursor,cols,where='',table='sdss',logx=False,logy=False):
    c,d=gettable(cursor,cols,where,table)
    if logx: c=N.log10(c)
    if logy: d=N.log10(d)
    P.plot(c,d,',',picker=5)
    P.title('X,Y: %s; Table: %s; Constr: %s'%(cols,table,where))

class inspect:
    def __init__(self,xcol,ycol,where='',curs=None,table='sdss',logx=False,logy=False,logxy=False,tolerance=0.01,smooth=4):
        if not curs: conn,curs=setupdb()
        self.xcol=xcol
        self.ycol=ycol
        #self.where=where+' AND %s NOT NULL AND %s NOT NULL'%(self.xcol,self.ycol)
        self.where=where
        self.curs=curs
        self.table=table
        self.smooth=smooth
        self.logx=logx
        self.logy=logy
        if logxy:
            self.logx=True
            self.logy=True
        self.tolerance=tolerance
        self.xtol=0.0
        self.ytol=0.0
        
        self.fig1=P.figure(1)
        self.init1()
        self.fig2=P.figure(2)
        self.init2()

        self.plot1()
        self.calctol()
        
        
    def init1(self):
        self.fig1.clf()
        self.ax1=self.fig1.add_subplot(1,1,1)
        self.ax1.set_xlabel(self.xcol)
        self.ax1.set_ylabel(self.ycol)
        self.ax1.set_title('Conditions: %s'%self.where)
        self.fig1.canvas.mpl_connect('key_press_event',self.keyhandler)
        self.clickconn=self.fig1.canvas.mpl_connect('pick_event',self.plotpicked)
        self.fig1.canvas.draw()


    def init2(self):
        self.fig2.clf()
        self.ax2=self.fig2.add_subplot(1,1,1)
        self.ax2.set_xlabel(u'Ångström')
        self.fig1.canvas.draw()

        
    def plot1(self):
        x,y=gettable(self.curs,self.xcol+','+self.ycol,self.where,self.table)
        if self.logx and self.logy: plotfu=self.ax1.loglog
        elif self.logx: plotfu=self.ax1.semilogx
        elif self.logy: plotfu=self.ax1.semilogy
        else: plotfu=self.ax1.plot
        plotfu(x,y,',',picker=5)
        self.fig1.canvas.draw()

    def calctol(self):
        axis=N.array(self.ax1.axis())
        self.xtol=self.tolerance*(axis[1]-axis[0])
        self.ytol=self.tolerance*(axis[3]-axis[2])
        
    def plotpicked(self,event):
        self.calctol()
        print self.xtol,self.ytol
        thisline = event.artist
        xdata,ydata = thisline.get_data()
        ind = event.ind
        x,y=zip(xdata[ind], ydata[ind])[0]
        where=self.where+' AND '
        where+='(%s BETWEEN %s AND %s) AND '%(self.xcol,x-self.xtol,x+self.xtol)
        where+='(%s BETWEEN %s AND %s)'%(self.ycol,y-self.ytol,y+self.ytol)
        ids=self.curs.execute('SELECT objID,z from %s WHERE %s'%(self.table,where))
        all=ids.fetchall()
        for id,z in all:
            #print id
            fits=specfromid(id,cursor=self.curs)
            head,spec,noise=sdss.splitfits(fits)
            wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
            wave/=1+z
            spec=smooth_gauss(spec,self.smooth)
            self.ax2.plot(wave,spec,linestyle='steps')

        self.fig2.canvas.draw()
        
            
    def keyhandler(self,event):
        if event.key=='s':
            print "Pick the point"
            self.clickconn=self.fig1.canvas.mpl_connect('pick_event',self.plotpicked)
        if event.key=='c': self.init2()
        if event.key=='q': self.fig1.canvas.mpl_disconnect(self.clickconn)
 
## specific plots below here

class inspectage(inspect):
    def __init__(self,xcol='age',ycol='Ha_w',where='z<5',curs=None,table='sb',logx=False,logy=False,logxy=True,tolerance=0.001,linef='/home/tom/projekte/sdss/ages/mixred0'):
        self.age,self.Ha_w=N.transpose(P.load(linef))[:2]
        inspect.__init__(self,xcol,ycol,where,curs,table,logx,logy,logxy,tolerance)
        
    def plot1(self):
        xcol=self.xcol.split(',')
        x1,y1=gettable(self.curs,xcol[0]+','+self.ycol,self.where,self.table)
        #x2,y2=gettable(self.curs,xcol[1]+','+self.ycol,self.where,self.table)
        #x3,y3=gettable(self.curs,xcol[2]+','+self.ycol,self.where,self.table)
        self.ax1.loglog(x1,y1,'g,',picker=5)
        #self.ax1.loglog(x2,y2,'r.',picker=5)
        #self.ax1.loglog(x3,y3,'b+',picker=5)
        self.ax1.loglog(self.age,self.Ha_w/2.0,'k--')
        self.ax1.loglog(self.age,self.Ha_w,'k-')
        self.fig1.canvas.draw()



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
    dage,dew,dmf,dM=gettable(cursor,cols='age,Ha_w,mf,M',where='agn=0 AND M>-20 AND age>5E5',table='sb')
    gage,gew,gmf,gM=gettable(cursor,cols='age,Ha_w,mf,M',where='agn=0 AND M<-20 AND age>1E6',table='sb')
    limit1x,limit1y=N.transpose(P.load('/home/tom/projekte/sdss/ages/mixred0'))[:2]
    limit2x,limit2y=N.transpose(P.load('/home/tom/projekte/sdss/ages/mixblue0'))[:2]
    #P.loglog(dage*2,dew,'g^',label='Age*2, M>-20',ms=5)
    #P.loglog(gage*2,gew,'b^',label='Age*2, M<-20',ms=5)
    
    P.subplot(1,2,1)
    P.loglog(limit1x,limit1y,'--r',linewidth=2)
    P.loglog(limit2x,limit2y,'-r',linewidth=2)
    P.scatter(dage*2,dew,c=dM,s=dmf*4000,label='Age*2, M>-20',alpha=0.2)
    P.xlabel(r'Burst age')
    P.ylabel(r'EW(H$\alpha$)')
    P.axis([1E6,2E10,100,4E3])
    P.legend()
    
    P.subplot(1,2,2)
    P.loglog(limit1x,limit1y,'--r',linewidth=2)
    P.loglog(limit2x,limit2y,'-r',linewidth=2)
    print N.log10(gmf)
    P.scatter(gage*2,gew,c=gM,s=gmf*4000,label='Age*2, M<-20',alpha=0.2)
    P.xlabel(r'Burst age')
    P.ylabel(r'EW(H$\alpha$)')
    P.axis([1E6,2E10,100,4E3])
    P.legend()

def plot2(cursor):
    extrawhere=' AND agn=0'
    plotbig,plotmid,plotsmall=True,True,True
    try: bigx,bigy=gettable(cursor,cols='NII_h/Ha_h,O5008_h/Hb_h',where='(Ha_s >700) AND (Hb_h >0) AND (O5008_h >0) AND (NII_h >0)'+extrawhere,table='sb')
    except: plotbig=False
    try: smallx,smally=gettable(cursor,cols='NII_h/Ha_h,O5008_h/Hb_h',where='(Ha_s BETWEEN 1 AND 250 ) AND (Hb_h >0) AND (O5008_h >0) AND (NII_h >0)'+extrawhere,table='sb')
    except: plotsmall=False
    try: midx,midy=gettable(cursor,cols='NII_h/Ha_h,O5008_h/Hb_h',where='(Ha_s BETWEEN 250 AND 700) AND (Hb_h >0) AND (O5008_h >0) AND (NII_h >0)'+extrawhere,table='sb')
    except: plotmid=False
    if plotbig: b=P.loglog(bigx,bigy,'^k',ms=5)
    if plotmid: m=P.loglog(midx,midy,'sg',ms=5)
    if plotsmall: s=P.loglog(smallx,smally,'Db',ms=5)
    #x=N.arange(0.001,0.4,0.01)
    #P.plot(x,10**sdss.stasinska(N.log10(x)),'r-',label='stasinska')
    x=N.arange(0.001,0.7,0.01)
    P.plot(x,10**sdss.mylee(N.log10(x)),'r-',linewidth=2)
    P.xlabel(r'$[N\, II]\, /\, H\alpha$')
    P.ylabel(r'$[O\, III]\, / H\beta$')
    #P.legend((s,m,b),(r'$\sigma(H\alpha)\, <\, 250\, km\,s^{-1}$',r'$250\, km\,s^{-1} <\, \sigma(H\alpha)\, <\, 700\, km\,s^{-1}$',r'$\sigma(H\alpha)\, >\, 700\, km\,s^{-1}$'),loc='lower left')

def plot3(curs):
    sbM,sbD,fade=gettable(curs,cols='M,voldens,fade',where='voldens NOTNULL AND M NOTNULL AND agn=0 AND age>5E5 AND bpara2 <3',table='sb')
    sbM2,sbD2,fade2=gettable(curs,cols='M,voldens,fade',where='voldens NOTNULL AND M NOTNULL AND agn=0 AND age>5E5 AND bpara2 >3',table='sb')
    pbM,pbD=gettable(curs,cols='M,voldens',where='voldens NOTNULL AND M NOTNULL',table='pb')
    abM,abD=gettable(curs,cols='M,voldens',where='voldens NOTNULL AND M NOTNULL AND agn=1 AND age>5E5',table='sb')
    sfM=N.array(sbM)+N.array(fade)
    X=N.arange(-24,-14,1.0,dtype='f')
    pby=sdss.lumfu(X,pbM,pbD)
    sby=sdss.lumfu(X,sbM,sbD)
    sby2=sdss.lumfu(X,sbM2,sbD2)
    sfy=sdss.lumfu(X,sfM,sbD)
    aby=sdss.lumfu(X,abM,abD)
    P.semilogy(X,sby2,'b-o',label=r'starbursts, EW(H$\alpha$)>120, b>3')
    P.semilogy(X,sby,'b--^',label=r'"starbursts", EW(H$\alpha$)>120, b<3')
    P.semilogy(X,pby,'r-o',label=r'postbursts, EW(H$\delta$)<-6)')
    P.semilogy(X[:5],aby[:5],'g-o',label='AGN')
    P.xlabel(r'$M_r$')
    P.ylabel(r'$\Phi$ [Mpc$^{-3}$]')
    P.legend(loc='lower right')

def plot4(curs):
    M,mgas,mstar,mtot=gettable(curs,cols='M,mgas,mass,mtot',where='mtot NOTNULL AND age > 5E5',table='sb')
    P.semilogy(M,mstar,'ro',label='mass in stars')
    P.semilogy(M,mtot*1E5,'go',label='total mass * 1E5')
    P.semilogy(M,mgas/1E5,'bo',label='mass in gas / 1E5')
    P.legend()

def plot5(curs):
    mtot,bpara,bpara2,age=gettable(curs,cols='mtot,bpara,bpara2,age',where='mtot NOTNULL AND age > 5E5 AND agn=0 AND bpara2 NOTNULL',table='sb')
    mtot=N.log10(mtot)
    X=N.arange(8,11,0.2,dtype='f')
    mean=sdss.averbins(X,mtot,bpara2)
    P.scatter(mtot,bpara2,s=age/1E8,label='b-parameter')
    P.plot(X[1:-1],mean[1:-1],'r-o')
    P.xlabel('total mass')
    P.ylabel('b-parameter')
    P.legend(('vertical mean','symbol size: age'))
    
def plot6(curs):
    sbM,sbD,fade,age=gettable(curs,cols='M,voldens,fade,age',where='voldens NOTNULL AND M NOTNULL AND agn=0 AND age>5E5 AND mf>0.01 AND age NOTNULL',table='sb')
    pbM,pbD=gettable(curs,cols='M,voldens',where='voldens NOTNULL AND M NOTNULL',table='pb')
    sfM=N.array(sbM)+N.array(fade)
    fact=8E8/age
    print fact
    sfD=sbD*fact
    X=N.arange(-24,-14,1.0,dtype='f')
    pby=sdss.lumfu(X,pbM,pbD)
    sby=sdss.lumfu(X,sbM,sbD)
    say=sdss.lumfu(X,sfM,sfD)
    sfy=sdss.lumfu(X,sfM,sbD)
    P.semilogy(X,pby,'r-o',label=r'postbursts')
    P.semilogy(X,sby,'b-o',label=r'starbursts, EW(H$\alpha$)>120, mf>1%')
    P.semilogy(X,sfy,'gD-',label=r'faded starbursts')
    P.semilogy(X,say,'r--^',label=r'faded starburts, age corrected')
    P.xlabel(r'$M_r$')
    P.ylabel(r'$\Phi$ [Mpc$^{-3}$]')
    P.legend(loc='lower right')

def plot7(curs):
    M,mgas,mstar,mtot,sfr=gettable(curs,cols='M,mgas,mass,mtot,sfr',where='mtot NOTNULL AND age>5E5 AND sfr NOTNULL',table='sb')
    P.loglog(mtot,mgas/sfr,'r,',label='mass in stars')
    P.xlabel(r'$M_{total}$')
    P.ylabel(r'$M_{gas} / SFR$')
    P.title('Gas consumption')

def plot8(curs):
    M1,mtot1,age1=gettable(curs,cols='M,mtot,age',where='mtot NOTNULL AND age>5E5 AND bpara2<3 AND mf <0.025 AND M NOTNULL',table='sb')
    M2,mtot2,age2=gettable(curs,cols='M,mtot,age',where='mtot NOTNULL AND age>5E5 AND bpara2<3 AND mf >0.025',table='sb')
    M3,mtot3,age3=gettable(curs,cols='M,mtot,age',where='mtot NOTNULL AND age>5E5 AND bpara2>3 AND mf <0.025',table='sb')
    M4,mtot4,age4=gettable(curs,cols='M,mtot,age',where='mtot NOTNULL AND age>5E5 AND bpara2>3 AND mf >0.025',table='sb')
    P.semilogy(M1,age1,'r.',label='b<3, mf <2.5%')
    P.semilogy(M2,age2,'k.',label='b<3, mf >2.5%')
    P.semilogy(M3,age3,'g.',label='b>3, mf <2.5%')
    P.semilogy(M4,age4,'b.',label='b>3, mf >2.5%')
    P.xlabel(r'$M_{i}$')
    P.ylabel(r'Age [yr]')
    P.legend(loc='lower right')

def plot9(curs):
    u,g,r = gettable(curs,cols='m_u,m_g,m_r',where='z<5',table='pb')

    P.plot(g-r,u-g,'ob')
    P.xlabel('g-r')
    P.ylabel('u-g')

def plot10(curs):
    Ha_h,Hb_h,z = gettable(curs,cols='Ha_h,Hb_h,z',where='z<5',table='sb')

    P.plot(z,Ha_h/Hb_h,'.b')
    P.grid()
    x=N.arange(0,0.4,0.001)
    y=sdss.averbins(x,z,Ha_h/Hb_h)
    P.plot(x,y,'or')
    P.xlabel('z')
    P.ylabel('Ha/Hb')
    

def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    

# TODO
# LumFu -> fading -> *alderskvot (8E8)
# gasmass/SFR vs log mass
# age burst vs absmag och massa, uppdelat i b>3 mf >2.5%.
