#!/usr/bin/env python
import numpy as np
import sqlite3
from pylab import *
from sdss import *
import pyfits
from itertools import izip

for t in (np.int8, np.int16, np.int32, np.int64,
        np.uint8, np.uint16, np.uint32, np.uint64):
        sqlite3.register_adapter(t, long)
for t in (np.float16, np.float32, np.float64, np.float128):
        sqlite3.register_adapter(t, double)

# data_mpa.sqlite should be fresh copy of data.sqlite
conn,curs = DB.setupdb('data_mpa.sqlite')

try:
    curs.execute("create table mpa (plate INTEGER, mjd INTEGER, fiber INTEGER, z FLOAT,EW_Ha FLOAT, EW_Hd FLOAT, balm_sigma FLOAT, SFR_fib FLOAT, SFR_tot FLOAT, m_stellar FLOAT);")
except:
    pass


ginfo=pyfits.open('mpa-sfr/gal_info_dr7_v5_2.fit.gz')[1].data
fibsfr=pyfits.open('mpa-sfr/gal_fibsfr_dr7_v5_2.fits.gz')[1].data
totsfr=pyfits.open('mpa-sfr/gal_totsfr_dr7_v5_2.fits.gz')[1].data
lines=pyfits.open('mpa-sfr/gal_line_dr7_v5_2.fit.gz')[1].data
stellar=pyfits.open('mpa-sfr/totlgm_dr7_v5_2.fit.gz')[1].data
#ourIDs=curs.execute('select plate,mjd,fiberID from sball').fetchall()

print "done opening files"

curs.executemany("insert into mpa values (?,?,?,?,?,?,?,?,?,?)",izip(\
            ginfo['PLATEID'],
            ginfo['MJD'],
            ginfo['FIBERID'],
            ginfo['Z'],
            lines['H_alpha_flux'] / lines['H_alpha_cont'],
            lines['H_delta_reqw']-lines['H_delta_eqw'],
            lines['sigma_balmer'],
            fibsfr['MEDIAN'],
            totsfr['MEDIAN'],
            stellar['MEDIAN']))
            #fibsfr['AVG'],
            #totsfr['AVG'],
            #stellar['AVG']))

curs.execute('create index platefib on mpa (plate,fiber);')
conn.commit()

conn.close()
