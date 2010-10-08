#!/usr/bin/env python

FIRSTFILE='objids.csv'

import sys
from sdss import *
from string import lower,strip

conn,curs=DB.setupdb()


try: f=open(sys.argv[1])
except:
    print "give the file to read!"
    exit()

if len(sys.argv) >2: SEP=sys.argv[2]
else: SEP=','
if SEP=='None':SEP=None

head=f.readline()
colist=map(strip,head.split(SEP))

if sys.argv[1]==FIRSTFILE:
    curs.execute('DROP TABLE IF EXISTS sdss;')
    curs.execute('CREATE TABLE sdss (objID INTEGER PRIMARY KEY, specObjID INTEGER);')
    curs.execute('create index objid_idx on sdss(objid);')
    curs.execute('create index specobjid_idx on sdss(specobjid);')
    qmarks=('?,'*len(colist))[:-1]
    sql="insert into sdss(%s) values (%s)"%(head,qmarks)
    print sql
    for line in f:
        values=map(strip,line.split(SEP))
        curs.execute(sql,values)


else:
    if lower(colist[0])=='specobjid': matchby='specObjID'
    else: matchby='objID'
    for col in colist: DB.createcolumnifnotexists(curs,col)
    setexp=('%s=?,'*len(colist))[:-1]
    setexp=setexp%tuple(colist)
    sql=" UPDATE sdss set %s WHERE %s=?"%(setexp,matchby)
    print sql
    for line in f:
        values=line[:-1].split(SEP)
        curs.execute(sql,values+[values[0]])



conn.commit()
conn.close()
    
