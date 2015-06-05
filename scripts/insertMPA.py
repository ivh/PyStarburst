#!/usr/bin/env python

from pylab import *
from sdss import *
import pyfits

conn,curs = DB.setupdb('data.sqlite')

#create table mpa (plate INTEGER, mjd INTEGER, fiber INTEGER, EW_Ha FLOAT, EW_Hd FLOAT, balm_sigma FLOAT, SFR_fib FLOAT, SFR_tot FLOAT);


ginfo=pyfits.open('mpa-sfr/gal_info_dr7_v5_2.fit.gz')[1].data
fibsfr=pyfits.open('mpa-sfr/gal_fibsfr_dr7_v5_2.fits.gz')[1].data
totsfr=pyfits.open('mpa-sfr/gal_totsfr_dr7_v5_2.fits.gz')[1].data
lines=pyfits.open('mpa-sfr/gal_line_dr7_v5_2.fit.gz')[1].data
ourIDs=curs.execute('select plate,mjd,fiberID from sball').fetchall()

print "done opening files"

count=0
for i,ginf in enumerate(ginfo):
    if tuple(ginf[:3]) in ourIDs:
        print ginf[:3], lines[i][:2]
        count+=1
        curs.execute("insert into mpa values (?,?,?,?,?,?,?,?)",(\
            int(ginf[0]),
            int(ginf[1]),
            int(ginf[2]),
            float(lines[i]['H_alpha_flux'] / lines[i]['H_alpha_cont']),
            float(lines[i]['H_delta_flux'] / lines[i]['H_delta_cont']),
            float(lines[i]['sigma_balmer']),
            float(fibsfr[i][0]),
            float(totsfr[i][0])))
        conn.commit()

print count
conn.close()
