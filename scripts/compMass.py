#!/usr/bin/env python

from pylab import *
from sdss import *
import pyfits

conn,curs = DB.setupdb('data_mpa.sqlite')


o_mass,t_mass,dyn1,dyn2,Mr,Mg,z,mz,mg,ML,MV=DB.get(curs,
    """select s.mtot, m.m_stellar,s.dynMassDisk,s.dynMassSphere,
     s.Mr,s.Mg,s.z,s.m_z,s.m_g,
     s.ML,s.MV_ec
     from sballz s, mpa m where s.fiberID=m.fiber and s.plate=m.plate and s.mjd=m.mjd""")



subplot(221)
hexbin(log10(o_mass),log10(dyn1),cmap=cm.bone_r)
X=arange(8,11,0.33)
#plot(sdss.averbins(X,log10(o_mass),log10(dyn1)),'ro-')
plot([0,700],[0,700],'r-')
axis([8,11,8,11])
xlabel('our Mtot')
ylabel('Mdyn')

subplot(222)
hexbin(t_mass,log10(dyn1),cmap=cm.bone_r)
X=arange(8,11,0.33)
#plot(sdss.averbins(X,log10(o_mass),log10(dyn1)),'ro-')
plot([0,700],[0,700],'r-')
axis([8,11,8,11])
xlabel('their Mtot')
ylabel('Mdyn')

subplot(223)
plot(MV,log10(ML),'k.',alpha=0.3)
axis([-14,-23,-1,1])
xlabel('M_V')
ylabel('our M/L')

print len(ML)

subplot(224)
ML_fac = (10**t_mass) / o_mass
t_ML=ML*ML_fac
t_ML=sdss.masked_where(t_ML<0,t_ML)
t_ML=sdss.masked_where(t_ML>5,t_ML)
plot(MV,log10(t_ML),'k.',alpha=0.3)
axis([-14,-23,-1,1])
xlabel('M_V')
ylabel('their M/L')



show()
conn.close()
