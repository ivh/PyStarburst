#!/usr/bin/env python

FIRSTFILE='objids.csv'

from sdss import *
import os
os.system('rm data.sqlite')
os.system('./csv2db.py objids.csv')
os.system('./csv2db.py spectradata2.csv')
os.system('./csv2db.py spectradata3.csv')
os.system('./csv2db.py photdata.csv')
os.system('./csv2db.py primtarget.csv')
os.system('./csv2db.py raddata.csv')
os.system('./csv2db.py sb_bestfits_20140325 None sbfit')
os.system('./csv2db.py best2z02_pb None pbfit')


conn,curs=DB.setupdb('data.sqlite')

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
S.fillDynMasses(curs)
S.fillBpara3(curs)

S.fillMass(curs,'pbfit','ML_tot') # postbursts as well
S.fillMtot(curs,'pbfit')

conn.commit()
conn.close()
