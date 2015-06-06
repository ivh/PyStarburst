#!/usr/bin/env python

from pylab import *
from sdss import *
import pyfits

conn,curs = DB.setupdb('data_mpa.sqlite')


o_Ha,t_Ha, t_sfr,o_sfr, o_mass,t_mass=DB.get(curs,
    'select s.EWHanew,m.EW_Ha,m.SFR_tot,s.SFR, m.m_stellar from sball s, mpa m where s.fiberID=m.fiber and s.plate=m.plate')

op_Hd,tp_Hd,op_mass,tp_mass = DB.get(curs,
    'select s.EW_Hd,m.EW_Hd, s.mtot,m.m_stellar from pball s, mpa m where s.fiberID=m.fiber and s.plate=m.plate')


subplot(221)
plot(o_Ha,t_Ha,'k.',alpha=0.2)
plot([50,700],[50,700],'r-')
axis([50,700,50,700])
title('H-aplha')

subplot(222)
plot(-op_Hd,tp_Hd,'k.',alpha=0.2)
plot([0,10],[0,10],'r-')
axis([5.5,11,0,5])
title('H-delta')

subplot(223)
plot(o_sfr,10**t_sfr,'k.',alpha=0.2)
plot([0,150],[0,150],'r-')
axis([0,150,0,150])
title('SFR')
xlabel('ours')
ylabel('theirs')

subplot(224)
plot(log10(o_mass),t_mass,'k.',alpha=0.2)
#plot(log10(op_mass),tp_mass,'r.',alpha=0.2)
plot([8,12],[8,12],'r-')
axis([8,12,8,12])
title('mass')


show()
conn.close()
