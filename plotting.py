#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import numpy as N
import pylab as P
from numpy.ma import masked_where

from db import *
import sdss
from PyGalKin.tool import smooth_gauss
from matplotlib.transforms import blended_transform_factory

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
        self.age,self.Ha_w=N.transpose(P.loadtxt(linef,unpack=True))[:2]
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



def plot1(cursor):
    dage,dew,dmf,dM=gettable(cursor,cols='age,Ha_w,mf,Mr',where='agn=0 AND Mr>-20 and mf NOTNULL',table='sb')
    gage,gew,gmf,gM=gettable(cursor,cols='age,Ha_w,mf,Mr',where='agn=0 AND Mr<-20 and mf NOTNULL',table='sb')
    limit1x,limit1y=P.loadtxt('/home/tom/projekte/sdss/ages/mixred0',unpack=True)[:2]
    limit2x,limit2y=P.loadtxt('/home/tom/projekte/sdss/ages/mixblue0',unpack=True)[:2]

    P.subplot(1,2,1)
    P.loglog(limit1x,limit1y,'--r',linewidth=2)
    P.loglog(limit2x,limit2y,'-r',linewidth=2)
    P.scatter(dage,dew,c=dM,s=dmf*2000,label='Age, M>-20',alpha=0.2)
    P.xlabel(r'Burst age')
    P.ylabel(r'EW(H$\alpha$)')
    P.axis([1.8E6,2E10,55,2E3])

    P.subplot(1,2,2)
    P.loglog(limit1x,limit1y,'--r',linewidth=2)
    P.loglog(limit2x,limit2y,'-r',linewidth=2)
    P.scatter(gage,gew,c=gM,s=gmf*2000,label='Age, M<-20',alpha=0.2)
    P.xlabel(r'Burst age')
    P.ylabel(r'EW(H$\alpha$)')
    P.axis([1.8E6,2E10,55,2E3])

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
    X=N.arange(-24,-15.5,1/3.0,dtype='f')
    P.semilogy(X,sdss.schechterBlanton(X),'k-',label='total (Blanton et al. (2001))')
    X=N.arange(-23,-15.5,1/3.0,dtype='f')

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND Ha_w > 100')
    y=sdss.lumfu(X,M,D)
    P.semilogy(X[2:],y[2:],'b-.s',label=r'$\mathrm{W(H\alpha)} > 100 \mathrm{\AA}$')

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND bpara2 >3')
    y=sdss.lumfu(X,M,D)
    P.semilogy(X,y,'b--*',label=r'$\mathrm{b} > 3$')

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND mf > 0.025')
    y=sdss.lumfu(X,M,D)
    P.semilogy(X,y,'b^:',label=r'mass fraction $> 2.5 \%$')

    M,D=getpb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL')
    y=sdss.lumfu(X,M,D)
    P.semilogy(X,y,'r-D',label=r'$\mathrm{W(H\delta)} < -6 \mathrm{\AA}$')

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=1')
    y=sdss.lumfu(X,M,D)
    y=masked_where(X>=-18,y)
    P.semilogy(X,y,'g-o',label='AGN')

    P.xlabel(r'$M_r$')
    P.ylabel(r'$\Phi\quad [\mathrm{Mpc}^{-3}\, \mathrm{mag}^{-1}]$')
    P.rcParams.update({'legend.fontsize':10})
    P.legend(loc='upper left')
    P.axis([-24.2,-15.5,4E-10,2E-2])

def plot4(curs):
    M,mgas,mstar,mtot=gettable(curs,cols='Mr,mgas,mass,mtot',where='mtot NOTNULL ',table='sb')
    P.semilogy(M,mstar,'ro',label='mass in stars')
    P.semilogy(M,mtot,'go',label='total mass')
    P.semilogy(M,mgas,'bo',label='mass in gas')
    P.legend()
    P.xlabel('Mr')

def plot5(curs):
    mtot,bpara,bpara2,age=gettable(curs,cols='mtot,bpara,bpara2,age',where='mtot NOTNULL AND agn=0 AND bpara NOTNULL',table='sb')
    b,label=bpara,'<b>'
    mtot=N.log10(mtot)
    X=N.arange(8,12,0.2,dtype='f')
    mean=sdss.averbins(X,mtot,b,median=True)
    P.scatter(mtot,b,s=age/1E8)
    P.plot(X[1:-1],mean[1:-1],'r-o')
    P.xlabel('log(total mass)')
    P.ylabel(label)
    P.legend(('median','symbol size: age'))

def plot6(curs):
    sbM,sbD,fade,age=gettable(curs,cols='Mr,voldens,fade,age',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND age>5E5 AND mf>0.05 AND age NOTNULL',table='sb')
    pbM,pbD=gettable(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL',table='pb')
    sfM=N.array(sbM)+N.array(fade)
    fact=8E8/age
    fusk=1
    sfD=sbD*fact
    X=N.arange(-24,-14,0.99,dtype='f')
    pby=sdss.lumfu(X,pbM,pbD)*fusk
    sby=sdss.lumfu(X,sbM,sbD)*fusk
    say=sdss.lumfu(X,sfM,sfD)*fusk
    sfy=sdss.lumfu(X,sfM,sbD)*fusk
    P.semilogy(X,pby,'r-o',label=r'postbursts')
    P.semilogy(X,sby,'b-o',label=r'starbursts, mf>5%')
    P.semilogy(X,sfy,'gD-',label=r'faded starbursts')
    P.semilogy(X,say,'r--^',label=r'faded starburts, age corrected')
    P.xlabel(r'$M_r$')
    P.ylabel(r'$\Phi$ [Mpc$^{-3}$]')
    P.legend(loc='lower right')
    P.axis([-24.5,-13.5,fusk*0.9E-13,fusk*1.05E-4])

def plot7(curs):
    M,mgas,mstar,mtot,sfr,bpara2=gettable(curs,cols='Mr,mgas,mass,mtot,sfr,bpara2',where='mtot NOTNULL AND sfr NOTNULL and agn=0',table='sb')
    P.loglog(mtot,mgas/sfr,'r.',label='b<3')
    mtot=masked_where(bpara2<3,mtot)
    P.loglog(mtot,mgas/sfr,'b.',label='b>3')
    P.xlabel(r'$M_{total}$')
    P.ylabel(r'$M_{gas} / SFR$')
    P.title('Gas consumption')
    P.legend(loc='upper left')


def plot8(curs):
    cols='Mr,mtot,age'
    where='mtot NOTNULL AND Mr NOTNULL'
    M1,mtot1,age1=gettable(curs,cols=cols,where=where+' AND mf <0.025 AND bpara<3',table='sb')
    M2,mtot2,age2=gettable(curs,cols=cols,where=where+' AND mf >0.025 AND bpara<3',table='sb')
    M3,mtot3,age3=gettable(curs,cols=cols,where=where+' AND mf <0.025 AND bpara>3',table='sb')
    M4,mtot4,age4=gettable(curs,cols=cols,where=where+' AND mf >0.025 AND bpara>3',table='sb')
    P.semilogy(M2,age2,'b.',label='b<3, mf >2.5%')
    P.semilogy(M1,age1,'r.',label='b<3, mf <2.5%')
#    P.semilogy(M3,age3,'g.',label='b>3, mf <2.5%')
#    P.semilogy(M4,age4,'b.',label='b>3, mf >2.5%')
    P.xlabel(r'$M_{r}$')
    P.ylabel(r'Age [yr]')
    P.legend(loc='lower right')

def plot8a(curs):
    cols='mtot,age'
    where='mtot NOTNULL AND Mr NOTNULL and bpara > 2.5'
    mtot1,age1=gettable(curs,cols=cols,where=where+' AND bpara2>3',table='sb')
    mtot2,age2=gettable(curs,cols=cols,where=where+' AND mf >0.025 AND bpara2>3',table='sb')
    mtot3,age3=gettable(curs,cols=cols,where=where+' AND mf >0.025',table='sb')
    P.subplot(131)
    P.ylabel(r'Age [yr]')
    P.title('b > 3')
    P.xlabel(r'total mass')
    P.loglog(mtot1,age1,'b.',label='b>3')
    axi=P.axis()

    P.subplot(132)
    P.title('b > 3, mf > 2.5%')
    P.loglog(mtot2,age2,'b.',label='b>3, mf >2.5%')
    ax=P.gca()
    ax.set_yticklabels([])
    P.axis(axi)

    P.subplot(133)
    P.title('mf >2.5%')
    P.loglog(mtot3,age3,'b.',label='mf >2.5%')
    ax=P.gca()
    ax.set_yticklabels([])
    P.axis(axi)

def plot9(curs):
    u,g,r = gettable(curs,cols='m_u,m_g,m_r',where='z<5',table='sb')

    P.plot(g-r,u-g,'ob')
    P.xlabel('g-r')
    P.ylabel('u-g')

def plot10(curs):
    Ha_h,Hb_h,z = gettable(curs,cols='Ha_h,Hb_h,z',where='Mr < -18 and Mr > -20',table='sb')

    P.plot(z,Ha_h/Hb_h,'.b')
    P.grid()
    x=N.arange(0,0.4,0.005)
    y=sdss.averbins(x,z,Ha_h/Hb_h)
    P.plot(x,y,'or')
    P.xlabel('z')
    P.ylabel('Ha/Hb')

def plot11(curs):
    age,mtot,mf,b=gettable(curs,cols='age,mtot,mf,bpara',where='bpara2 > 3',table='sb')
    P.loglog(mtot,age,'Db',label='b>3')
    mtot=masked_where(mf<0.025,mtot)
    P.loglog(mtot,age,'Dg',label='b>3, mf>2.5%')
    mtot=masked_where(b<3,mtot)
    P.loglog(mtot,age,'Dr',label='b>3, mf>2.5%, <b> >3')
    P.xlabel('total mass')
    P.ylabel('age')
    P.legend(loc='lower right')

def plot12(curs):
    age,mtot,b=gettable(curs,cols='age,mtot,bpara',where='agn=0',table='sb')
    b1=masked_where(age>=5E7,b)
    b2=masked_where(age<=5E7,b)
    b2=masked_where(age>=5E8,b2)
    b3=masked_where(age<=5E8,b)
    P.semilogx(mtot,b3,',r',label='age > 5E8 yr')
    P.semilogx(mtot,b2,'.b',label='5E7 < age < 5E8 yr')
    P.semilogx(mtot,b1,'.k',label='age < 5E7 yr')
    P.axis([5E8,1E12,-1,15])
    P.xlabel('total mass')
    P.ylabel('<b>')
    P.legend(loc='upper right')

def plot13(curs):
    b,Ha_w,Ha_w_fit=gettable(curs,cols='bpara2,Ha_w,Ha_w_fit',where='Ha_w > 120',table='sb')
    P.plot(Ha_w_fit,b,'r.',label='Ha_w_fit')
    P.plot(Ha_w,b,'b.',label='Ha_w')
    P.xlabel('W(Ha)')
    P.ylabel('b')
    P.axis([0,500,-1,10])
    P.legend(loc='upper right')

def plot14(curs):
    dynMassDisk,dynMassSphere=gettable(curs,cols='dynMassDisk,dynMassSphere',where='agn=0',table='sb')
    P.loglog(dynMassDisk,dynMassSphere,'.b')
    dynMassDisk,dynMassSphere=gettable(curs,cols='dynMassDisk,dynMassSphere',where='agn=0',table='pb')
    P.loglog(dynMassDisk,dynMassSphere,'.r')
    P.grid()
    P.xlabel('M_dyn_disk')
    P.ylabel('M_dyn_sphere')

def plot15(curs):
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0',table='sb')
    md=(md1+md2)/2.0
    P.loglog(mp,md,'.b')
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0',table='pb')
    md=(md1+md2)/2.0
    P.loglog(mp,md,'.r')
    P.grid()
    P.xlabel('M_phot')
    P.ylabel('M_dyn')

def plot16(curs):
    b1,b2=gettable(curs,cols='bpara2,bpara3',where='agn=0',table='sb')
    P.plot(b1,b2,'.b')
    b1,b2=gettable(curs,cols='bpara2,bpara3',where='agn=0',table='pb')
    P.plot(b1,b2,'.r')
    P.grid()
    P.xlabel('b_phot')
    P.ylabel('b_dyn')
    P.axis([-1,20,-1,20])

def plot17(curs):
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0 and mtot NOTNULL',table='sb')
    md=(md1+md2)/2.0
    P.semilogx(mp,md/mp,'.b')
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0 and mtot NOTNULL',table='pb')
    md=(md1+md2)/2.0
    P.semilogx(mp,md/mp,'.r')
    P.grid()
    P.xlabel('M_phot')
    P.ylabel('M_dyn/M_phot')

def plot18(curs):
    bins=N.array([5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0])
    age1=getsb(curs,cols='age',where='agn=0 and Ha_w > 100 and age NOTNULL')
    age2=getsb(curs,cols='age',where='agn=0 and bpara2 > 3 and age NOTNULL')
    age3=getsb(curs,cols='age',where='agn=0 and mf > 0.025 and age NOTNULL')
    ages=map(N.log10,(age1,age2,age3))
    titles=[r'$\mathrm{W(H\alpha)}\, >\, 100 \mathrm{\AA}$',
            r'$\mathrm{b}\, >\, 3$',
            r'$\mathrm{mass\, fraction}\, >\, 2.5 \%$']
    for i,age in enumerate(ages):
        ax=P.subplot(131+i)
        P.hist(x=age,normed=True,bins=bins,fc='gray')
        P.title(titles[i])
        P.xlabel(r'$\mathrm{log}_{10}(\mathrm{burst\,  age})$')
        P.xticks(N.arange(4)+6)
        P.yticks([])
        trans=blended_transform_factory(ax.transData, ax.transAxes)
        P.plot([7.0],[0.8],'k<', transform=trans)
        P.plot([8.0],[0.8],'ko', transform=trans)
        P.plot([8.65],[0.8],'k>', transform=trans)

def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
