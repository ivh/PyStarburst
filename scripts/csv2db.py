#!/usr/bin/env python

FIRSTFILE='objids.csv'

import sys
from sdss import *
from string import lower,strip

conn,curs=DB.setupdb('data.sqlite')


try: f=open(sys.argv[1])
except:
    print "give the file to read!"
    exit()

if len(sys.argv) >2: SEP=sys.argv[2]
else: SEP=','
if SEP=='None':SEP=None

if len(sys.argv) >3: TABLE = sys.argv[3]
else: TABLE='sdss'

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


elif TABLE=='sdss':
    if lower(colist[0])=='specobjid': matchby='specObjID'
    else: matchby='objID'
    colist = colist[1:] # REMOVE objid or specobjid
    for col in colist: DB.createcolumnifnotexists(curs,col)
    setexp=('%s=?,'*len(colist))[:-1]
    setexp=setexp%tuple(colist)
    sql=" UPDATE sdss set %s WHERE %s=?"%(setexp,matchby)
    print sql
    for line in f:
        values=line.strip().split(SEP)
        curs.execute(sql,values[1:]+[values[0]])

else: # we read into a new table
    print 'making new table %s'%TABLE
    firstline=f.readline()
    coldef=' '
    for i,val in enumerate(firstline.split(SEP)[1:]):
        try:
            f=float(val)
            if float(int(val))==f:
                coldef+='%s INTEGER, '%colist[i+1]
            else: coldef+='%s REAL, '%colist[i+1]
        except:
            coldef+='%s TEXT, '%colist[i+1]

    print colist,coldef
    curs.execute('CREATE TABLE IF NOT EXISTS %s (%s INTEGER PRIMARY KEY, %s)'%\
        (TABLE,colist[0],coldef[:-2]))

    qmarks=('?,'*len(colist))[:-1]
    sql="insert into %s values (%s)"%(TABLE,qmarks)
    values=map(strip,firstline.split(SEP))
    curs.execute(sql,values)
    for line in f:
        values=map(strip,line.split(SEP))
        curs.execute(sql,values)


conn.commit()
conn.close()
