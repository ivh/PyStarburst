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
#from PyGalKin.tool import smooth_gauss
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

def plot_zhist(curs):
    sbz,sbb=gettable(curs,cols='z,b_para',where='agn=0',table='sball')
    pbz,pbb=gettable(curs,cols='z,b_param',where='z>0.0',table='pball')
    bins=N.arange(0.0,0.24,0.02)
    P.hist(sbz,bins=bins,ls='solid',color='k',histtype='step',label=r'starburst candidates, $\mathrm{EW(H_\alpha )}>60\mathrm{\AA}$')
    P.hist(pbz,bins=bins,ls='dashed',color='r',histtype='step',label='postburst candidates, $\mathrm{EW(H_\delta )}<-6\mathrm{\AA}$')

    sbz,sbb=gettable(curs,cols='z,b_para',where='agn=0 and Massfrac > 0.03',table='sball')
    pbz,pbb=gettable(curs,cols='z,b_param',where='Massfrac > 0.03',table='pball')
    P.hist(pbz,bins=bins,histtype='stepfilled',color='r',alpha=0.3,label=r'postbursts $f_{burst}>3\%$',weights=N.ones(pbz.shape))
    P.hist(sbz,bins=bins,histtype='stepfilled',color='k',alpha=0.3,label=r'starbursts $f_{burst}>3\%$',weights=N.ones(sbz.shape))


    P.legend(loc='upper right',fontsize=11)
    P.xlabel(r'$z$')
    P.ylabel(r'$N$')


def plot1(cursor):
    dage,dew,dmf,dM=gettable(cursor,cols='age,Ha_w,Massfrac,Mr',where='agn=0 AND Mr>-20 and Massfrac NOTNULL',table='sball')
    gage,gew,gmf,gM=gettable(cursor,cols='age,Ha_w,Massfrac,Mr',where='agn=0 AND Mr<-20 and Massfrac NOTNULL',table='sball')
    limit1x,limit1y=P.loadtxt('/home/tom/projekte/sdss/ages/mixred0',unpack=True)[:2]
    limit2x,limit2y=P.loadtxt('/home/tom/projekte/sdss/ages/mixblue0',unpack=True)[:2]

    ax1=P.subplot(1,2,1)
    P.title(r'$\mathrm{Dwarfs,}\,\, M_r > -20$')
    P.loglog(limit1x,limit1y,'--k',linewidth=2)
    P.loglog(limit2x,limit2y,'-k',linewidth=2)
    P.scatter(dage,dew*dcf,c=dmf,vmin=0,vmax=0.05,label='Age, M>-20',alpha=0.2,cmap=P.cm.OrRd,lw=1)
    P.xlabel(r'$\mathrm{Burst\, age}$')
    P.ylabel(r'$\mathrm{EW}(H\alpha)$')
    P.plot([1E6,1E10],[150,150],'k-.',lw=2)
    P.axis([4E6,2E9,55,9E2])

    P.subplot(1,2,2)
    P.title(r'$M_r < -20$')
    P.loglog(limit1x,limit1y,'--k',linewidth=2)
    P.loglog(limit2x,limit2y,'-k',linewidth=2)
    P.scatter(gage,gew*gcf,c=gmf,vmin=0,vmax=0.05,label='Age, M<-20',alpha=0.2,cmap=P.cm.OrRd,lw=1)
    P.xlabel(r'$\mathrm{Burst\, age}$')
    P.plot([1E6,1E10],[150,150],'k-.',lw=2)
    P.axis([4E6,2E9,55,9E2])
    P.setp(P.gca(),'yticklabels',[])
    P.subplots_adjust(wspace=0)

def plot2(cursor):
    ax,ay=gettable(cursor,cols='NII_h/Ha_h,O5008_h/Hb_h',where='agn=1',table='sb')
    nx,ny=gettable(cursor,cols='NII_h/Ha_h,O5008_h/Hb_h',where='agn=0',table='sb')
    P.loglog(ax,ay,',g',ms=5)
    P.loglog(nx,ny,'Db',ms=5,alpha=0.5)
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

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND Massfrac> 0.03')
    y=sdss.lumfu(X,M,D)
    P.semilogy(X,y,'b^:',label=r'mass fraction $> 3 \%$')

    M,D=getpb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL')
    y=sdss.lumfu(X,M,D)
    P.semilogy(X,y,'r-D',label=r'$\mathrm{W(H\delta)} < -6 \mathrm{\AA}$')

    M,D=gettable(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=1',table='sb')
    y=sdss.lumfu(X,M,D)
    y=masked_where(X>=-18,y)
    P.semilogy(X,y,'g-o',label='AGN')

    P.xlabel(r'$M_r$')
    P.ylabel(r'$\Phi\quad [\mathrm{Mpc}^{-3}\, \mathrm{mag}^{-1}]$')
    P.rcParams.update({'legend.fontsize':10})
    P.legend(loc='upper left')
    P.axis([-24.2,-15.5,4E-10,2E-2])

def plot4(curs):
    M,mgas,mstar,mtot,msph,mdisk=gettable(curs,cols='Mr,mgas,mass,mtot,dynMassDisk,dynMassSphere',where='mtot NOTNULL ',table='sball')
    P.semilogy(M,msph,'yo',label='dyn. mass (sphere)')
    P.semilogy(M,mdisk,'ko',label='dyn. mass (disk)')
    P.semilogy(M,mstar,'ro',label='mass in stars')
    P.semilogy(M,mtot,'go',label='total mass')
    P.semilogy(M,mgas,'bo',label='mass in gas')
    P.legend()
    P.xlabel(r'$M_r$')
    P.ylabel(r'$Mass\quad(M_\odot)$')
    P.title(r'Mass comparisons (starbursts only)')

def plot5(curs):
    f=P.gcf()
    ax1=f.add_axes((.12,.53,.85,.45))
    ax2=f.add_axes((.12,.09,.85,.45))
    mtot,massfrac=gettable(curs,cols='mtot,massfrac',where='mtot NOTNULL AND agn=0 AND massfrac NOTNULL AND bpara2 >3',table='sball')
    mtot=N.log10(mtot)
    X=N.arange(8,11.4,0.2,dtype='f')
    mean=sdss.averbins(X,mtot,massfrac,median=True)
    ax1.hexbin(mtot,massfrac,cmap=P.cm.bone_r,gridsize=(50,80))
    ax1.plot(X[4:-1],mean[4:-1],'k-o',lw=2, markeredgecolor='w', markeredgewidth=1.5)
    #ax1.plot(X[4:-1],mean[4:-1],'',lw=12)

    #massfrac=masked_where(massfrac<0.03,massfrac)
    #mean=sdss.averbins(X,mtot,massfrac,median=True)
    #ax1.plot(X[4:-1],mean[4:-1],'y--D')

    ax1.set_ylabel(r'$\mathrm{mass\,\,fraction}$')
    ax1.axis((8.4,11.2,-.001,0.15))
    ax1.text(8.7,0.13,r'$\mathrm{starbursts,\,\, b > 3}$',fontsize=11)

    mtot,massfrac=gettable(curs,cols='mtot,massfrac',where='mtot NOTNULL AND agn=0 AND massfrac NOTNULL',table='pball')
    mtot=N.log10(mtot)
    X=N.arange(8,11.4,0.2,dtype='f')
    mean=sdss.averbins(X,mtot,massfrac,median=True)
    ax2.hexbin(mtot,massfrac,cmap=P.cm.bone_r,gridsize=(50,80))
    ax2.plot(X[4:-1],mean[4:-1],'k-o',lw=2, markeredgecolor='w', markeredgewidth=1.5)
    ax2.set_ylabel(r'$\mathrm{mass\,\,fraction}$')
    ax2.set_xlabel(r'$\log_{10}({\cal M}_{tot})$')
    ax2.axis((8.4,11.2,-.001,0.15))
    ax2.text(8.7,0.13,r'$\mathrm{postbursts}$',fontsize=11)

def plot6(curs):
    sbM,sbD,fade,age=gettable(curs,cols='Mr,voldens,fade,age',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND Massfrac>0.025 AND age NOTNULL',table='sball')
    pbM,pbD=gettable(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL',table='pball')
    sfM=N.array(sbM)+N.array(fade)
    fact=8E8/age
    fusk=1
    sfD=sbD*fact
    X=N.arange(-23,-15.5,1/3.0,dtype='f')
    pby=sdss.lumfu(X,pbM,pbD)*fusk
    sby=sdss.lumfu(X,sbM,sbD)*fusk
    say=sdss.lumfu(X,sfM,sfD)*fusk
    sfy=sdss.lumfu(X,sfM,sbD)*fusk
    P.semilogy(X,pby,'r-o',label=r'postbursts')
    P.semilogy(X,sby,'b-o',label=r'starbursts, $mf>2.5\%$')
    P.semilogy(X,sfy,'gD-',label=r'faded starbursts')
    P.semilogy(X,say,'r--^',label=r'faded starburts, age corrected')
    P.xlabel(r'$M_r$')
    P.ylabel(r'$\Phi\quad [\mathrm{Mpc}^{-3}\, \mathrm{mag}^{-1}]$')
    P.rcParams.update({'legend.fontsize':10})
    P.legend(loc='upper left')
    P.axis([-24.2,-15.5,4E-10,2E-2])

def plot7(curs):
    M,mgas,mstar,mtot,sfr,bpara2=gettable(curs,cols='Mr,mgas,mass,mtot,sfr,bpara2',where='mtot NOTNULL AND sfr NOTNULL and agn=0',table='sball')
    P.loglog(mtot,mgas/sfr,'r.',label=r'$b<3$',alpha=0.3)
    mtot=masked_where(bpara2<3,mtot)
    P.loglog(mtot,mgas/sfr,'b.',label=r'$b>3$',alpha=0.3)
    P.xlabel(r'$M_{total}$')
    P.ylabel(r'$M_{gas} / SFR$')
    P.title('Gas consumption')
    P.legend(loc='upper left')
    P.axis([1E7,1E12,1E8,1E11])

def plot8(curs):
    cols='mtot,age'
    where='mtot NOTNULL AND Mr NOTNULL and b_para > 2.5'
    mtot1,age1=gettable(curs,cols=cols,where=where+' AND bpara2>3',table='sball')
    mtot2,age2=gettable(curs,cols=cols,where=where+' AND Massfrac >0.025 AND bpara2>3',table='sball')
    mtot3,age3=gettable(curs,cols=cols,where=where+' AND Massfrac >0.025',table='sball')
    P.subplot(131)
    P.ylabel(r'Age [yr]')
    P.title(r'$b > 3$')
    P.xlabel(r'$m_{tot}$')
    P.loglog(mtot1,age1,'b.',label='b>3')
    axi=P.axis()

    P.subplot(132)
    P.title(r'$b > 3, mf > 2.5\%$')
    P.loglog(mtot2,age2,'b.',label='b>3, mf >2.5%')
    ax=P.gca()
    ax.set_yticklabels([])
    P.axis(axi)

    P.subplot(133)
    P.title(r'$mf >2.5\%$')
    P.loglog(mtot3,age3,'b.',label='mf >2.5%')
    ax=P.gca()
    ax.set_yticklabels([])
    P.axis(axi)

def plot9(curs):
    ax1=P.subplot(1,2,1)
    u,g,r = gettable(curs,cols='m_u,m_g,m_r',where='agn=1',table='sb')
    P.plot(g-r,u-g,'.g',alpha=0.2,label='agn')
    u,g,r = gettable(curs,cols='m_u,m_g,m_r',where='agn=0',table='sb')
    P.plot(g-r,u-g,'.b',alpha=0.3,label='starbursts')
    P.axis([-1,2,-0.5,3])
    P.xlabel('g-r')
    P.ylabel('u-g')
    ax2=P.subplot(1,2,2)
    u,g,r = gettable(curs,cols='m_u,m_g,m_r',where='agn=1',table='pb')
    P.plot(g-r,u-g,'.r',alpha=0.3,label='postbursts')
    P.axis([-1,2,-0.5,3])
    P.setp(P.gca(),'yticklabels',[])
    P.subplots_adjust(wspace=0)


def plot10(curs):
    Ha_h,Hb_h,z = gettable(curs,cols='Ha_h,Hb_h,z',where='Mr < -18 and Mr > -20',table='sball')

    P.plot(z,Ha_h/Hb_h,'.b')
    P.grid()
    x=N.arange(0,0.4,0.005)
    y=sdss.averbins(x,z,Ha_h/Hb_h)
    P.plot(x,y,'or')
    P.xlabel('z')
    P.ylabel('Ha/Hb')

def plot11(curs):
    age,mtot,mf,b=gettable(curs,cols='age,mtot,Massfrac,b_para',where='bpara2 > 3',table='sball')
    P.loglog(mtot,age,'Db',label=r'$b>3$')
    mtot=masked_where(mf<0.025,mtot)
    P.loglog(mtot,age,'Dg',label=r'$b>3, mf>2.5\%$')
    mtot=masked_where(b<3,mtot)
    P.loglog(mtot,age,'Dr',label=r'$b>3, mf>2.5\%, <b> >3$')
    P.xlabel(r'$m_{tot}$')
    P.ylabel('age')
    P.legend(loc='lower right')

def plot12(curs):
    age,mtot,b=gettable(curs,cols='age,mtot,b_para',where='agn=0',table='sball')
    b1=masked_where(age>=5E7,b)
    b2=masked_where(age<=5E7,b)
    b2=masked_where(age>=5E8,b2)
    b3=masked_where(age<=5E8,b)
    P.semilogx(mtot,b3,',r',label=r'age $> 5E8 yr$')
    P.semilogx(mtot,b2,'.b',label=r'$5E7 <$ age $< 5E8 yr$')
    P.semilogx(mtot,b1,'.k',label=r'age $< 5E7 yr$')
    P.axis([5E8,1E12,-1,15])
    P.xlabel(r'$m_{tot}$')
    P.ylabel(r'$<b>$')
    P.legend(loc='upper right')

def plot13(curs):
    b,Ha_w,Ha_w_fit,age=gettable(curs,cols='bpara2,Ha_w,EWHanew,age',where='Ha_w > 120',table='sball')
    P.plot(Ha_w_fit,N.ma.masked_where(age < 1E8,b),'r.',label=r'$\mathrm{age} > 1E8$', alpha=0.5)
    P.plot(Ha_w_fit,N.ma.masked_where(age > 1E8,b),'b.',label=r'$\mathrm{age} < 1E8$', alpha=0.5)
    #P.plot(Ha_w,b,'b.',label=r'$\mathrm{EW}(H\alpha)$', alpha=0.5)
    P.xlabel(r'$\mathrm{EW}(H\alpha)$')
    P.ylabel(r'$b$')
    P.axis([0,500,-1,10])
    P.legend(loc='upper right')

def plot14(curs):
    dynMassDisk,dynMassSphere=gettable(curs,cols='dynMassDisk,dynMassSphere',where='agn=0',table='sball')
    P.loglog(dynMassDisk,dynMassSphere,'.b', alpha=0.5, label='starbursts')
    dynMassDisk,dynMassSphere=gettable(curs,cols='dynMassDisk,dynMassSphere',where='agn=0',table='pb')
    P.loglog(dynMassDisk,dynMassSphere,'.r', alpha=0.5, label='postbursts')
    P.grid()
    P.xlabel(r'$M_{dyn,disk}$')
    P.ylabel(r'$M_{dyn,sphere}$')

def plot15(curs):
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0',table='sball')
    md=(md1+md2)/2.0
    P.loglog(mp,md,'.b',alpha=0.5)
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0',table='pball')
    md=(md1+md2)/2.0
    P.loglog(mp,md,'.r')
    P.grid()
    P.xlabel(r'$M_{phot}$')
    P.ylabel(r'$M_{dyn}$')

def plot16(curs):
    b1,b2=gettable(curs,cols='bpara2,bpara3',where='agn=0',table='sball')
    P.plot(b1,b2,'.b',alpha=0.5)
    #b1,b2=gettable(curs,cols='bpara2,bpara3',where='agn=0',table='pball')
    #P.plot(b1,b2,'.r')
    P.grid()
    P.xlabel(r'$b_{phot}$')
    P.ylabel(r'$b_{dyn}$')
    P.axis([-0.1,5,-0.1,5])

def plot17(curs):
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0 and mtot NOTNULL',table='sball')
    md=(md1+md2)/2.0
    P.semilogx(mp,md/mp,'.b',alpha=0.5)
    mp,md1,md2=gettable(curs,cols='mtot,dynMassDisk,dynMassSphere',where='agn=0 and mtot NOTNULL',table='pball')
    md=(md1+md2)/2.0
    P.semilogx(mp,md/mp,'.r')
    P.grid()
    P.xlabel(r'$M_{phot}$')
    P.ylabel(r'$M_{dyn} / M_{phot}$')
    P.axis([2E7,9E11,-1,80])

def plot18(curs):
    bins=N.array([5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0])
    age1=getsb(curs,cols='age',where='agn=0 and Ha_w > 100 and age NOTNULL')
    age2=getsb(curs,cols='age',where='agn=0 and bpara2 > 3 and age NOTNULL')
    age3=getsb(curs,cols='age',where='agn=0 and Massfrac > 0.025 and age NOTNULL')
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

def plot19(curs):
    Fig=P.gcf()
    ax1=Fig.add_axes((.12,.4,.85,.58))
    ax2=Fig.add_axes((.12,.07,.85,.33))
    X=N.arange(-24,-15.5,1/3.0,dtype='f')
    Blanton=sdss.schechterBlanton(X)
    ax1.semilogy(X,Blanton,'k-',label='total (Blanton et al. (2001))')

    # completeness limit
    ax1.plot((-17,-17),(1E-8,1E-3),'k--',lw=3)
    ax2.plot((-17,-17),(1E-4,2E-2),'k--',lw=3)

    #M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND Ha_w > 100')
    #y=sdss.lumfu(X,M,D)
    #ax1.semilogy(X,y,'b-.s',label=r'$\mathrm{W(H\alpha)} > 100 \mathrm{\AA}$')
    #ax2.semilogy(X,y/Blanton,'b-.s',label=r'$\mathrm{W(H\alpha)} > 100 \mathrm{\AA}$')

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND bpara2 >3')
    y=sdss.lumfu(X,M,D)
    ax1.semilogy(X,y,'b--*',label=r'$\mathrm{b} > 3$')
    ax2.semilogy(X,y/Blanton,'b--*',label=r'$\mathrm{b} > 3$')

    M,D=getsb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=0 AND Massfrac> 0.03')
    y=sdss.lumfu(X,M,D)
    ax1.semilogy(X,y,'b^:',label=r'mass fraction $> 3 \%$')
    ax2.semilogy(X,y/Blanton,'b^:',label=r'mass fraction $> 3 \%$')

    M,D=getpb(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL')
    y=sdss.lumfu(X,M,D)
    ax1.semilogy(X,y,'r-D',label=r'$\mathrm{W(H\delta)} < -6 \mathrm{\AA}$')
    ax2.semilogy(X,y/Blanton,'r-D',label=r'$\mathrm{W(H\delta)} < -6 \mathrm{\AA}$')

    M,D=gettable(curs,cols='Mr,voldens',where='voldens NOTNULL AND Mr NOTNULL AND agn=1',table='sb')
    y=sdss.lumfu(X,M,D)
    y=masked_where(X>=-18,y)
    ax1.semilogy(X[3:],y[3:],'g-o',label='AGN')
    ax2.semilogy(X[3:],(y/Blanton)[3:],'g-o',label='AGN')

    ax1.set_xticklabels([])
    ax1.set_ylabel(r'$\Phi\quad [\mathrm{Mpc}^{-3}\, \mathrm{mag}^{-1}]$')
    P.rcParams.update({'legend.fontsize':10})
    ax1.legend(loc='upper left')
    ax1.axis([-24.2,-15.5,4E-10,2E-2])

    ax2.set_xlabel(r'$M_r$')
    ax2.set_ylabel(r'$\mathrm{Fraction}$')
    ax2.axis([-24.2,-15.5,4E-5,4E-2])


def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
