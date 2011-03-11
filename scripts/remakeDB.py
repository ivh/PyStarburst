#!/usr/bin/env python

FIRSTFILE='objids.csv'

from sdss import *
import os
os.system('rm sdss.db')
os.system('./csv2db.py objids.csv')
os.system('./csv2db.py spectradata2.csv')
os.system('./csv2db.py spectradata3.csv')
os.system('./csv2db.py photdata.csv')
os.system('./csv2db.py primtarget.csv')
os.system('./csv2db.py raddata.csv')
os.system('./csv2db.py sbresult None')
# old data
#os.system('./csv2db.py mix100901 None')


conn,curs=DB.setupdb()

S.createviews(curs)
conn.commit()

S.fixSigmaInDB(curs) # Angstrom to Km/s
conn.commit()

S.decideAGN(curs)
S.fillM(curs)
S.fillVoldens(curs)
S.fillMass(curs)
S.fillFading(curs)
S.fillSFR(curs)
S.fillMgas(curs)
S.fillMtot(curs)
S.fillBpara2(curs)


conn.commit()
conn.close()
    
