#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import numpy as N
import pylab as P

from db import *
from sdss import splitfits
from PyGalKin.tool import smooth_gauss


def plotspecwhere(cursor,where,table='sdss'):
    ids=cursor.execute('SELECT objID from %s WHERE %s'%(table,where))
    for id in ids:
       plotspecbyid(id) 

       
def plotspecbyid(id):
    fits=specfromid(id)
    head,spec,noise=splitfits(fits)
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
            head,spec,noise=splitfits(fits)
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
    def __init__(self,xcol='age3',ycol='Ha_w',where='z<5',curs=None,table='sb',logx=False,logy=False,logxy=True,tolerance=0.001,linef='/home/tom/projekte/sdss/ages/mixred0'):
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
    big=getsb(cursor,'Ha_s,Ha_h,NII_h','Ha_s >15')
    small=getsb(cursor,'Ha_s,Ha_h,NII_h','Ha_s <5')
    midd=getsb(cursor,'Ha_s,Ha_h,NII_h','(Ha_s >5 AND Ha_s <15)')
    numsb=midd[0].size+big[0].size+small[0].size
    P.plot(N.log10(big[2]/big[1]),N.log10(big[0]),',k')
    P.plot(N.log10(small[2]/small[1]),N.log10(small[0]),',b')
    P.plot(N.log10(midd[2]/midd[1]),N.log10(midd[0]),',g')
    P.title('Starbursts (%d)'%numsb)
    P.xlabel('log NII/Ha')
    P.ylabel('log sigma Ha')

def plot2(cursor):
    big=getsb(cursor,'O5008_h/Hb_h,Ha_h,NII_h','Ha_s >15')
    small=getsb(cursor,'O5008_h/Hb_h,Ha_h,NII_h','Ha_s <5')
    midd=getsb(cursor,'O5008_h/Hb_h,Ha_h,NII_h','(Ha_s >5 AND Ha_s <15)')
    numsb=midd[0].size+big[0].size+small[0].size
    P.plot(N.log10(big[2]/big[1]),N.log10(big[0]),',k')
    P.plot(N.log10(small[2]/small[1]),N.log10(small[0]),',b')
    P.plot(N.log10(midd[2]/midd[1]),N.log10(midd[0]),',g')
    P.title('Starbursts, %s galaxies'%numsb)
    P.xlabel('log NII/Ha')
    P.ylabel('log O5008 / Hb')

    
def demo():
    print "This file defines some functions. It is not meant to be executed. Import it instead!"

if __name__ == '__main__':
    demo()
    
