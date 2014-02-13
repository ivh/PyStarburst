#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help:
* call the script with all the csv-files that you
  want to have in your DB.
* first line has to contain the names of the columns
* second line has to contain the data type of each
  column, e.g. INTEGER, REAL, TEXT.
* each file has to have a column called objID which
  is used to identify the rows.
* the database is created if not existing. the script
  inserts into existing databases.
* change the path to the database file on top of the code
  or give a filename ending with .db as first argument
"""

# IMPORTS
import sys
import os
import string
from os.path import join,exists
from sqlite3 import dbapi2 as sqlite
import numpy as N
import pylab as P
import pyfits as F
import sdss

# GLOBAL VARIABLES
DBNAME=join('/','data','sdss','sdss.db')
DELIM=','
SPECBASE=join('/','data','sdss','spectra')
SPECZBASE=join('/','data','sdss','spZline')

# The list of lines in the SpZ-fitsfiles
SpZlist=['Ly_alpha','N_V 1240','C_IV 1549','He_II 1640','C_III] 1908','Mg_II 2799','[O_II] 3725','[O_II] 3727','[Ne_III] 3868','H_epsilon','[Ne_III] 3970','H_delta','H_gamma','[O_III] 4363','He_II 4685','H_beta','[O_III] 4959','[O_III] 5007','He_II 5411','[O_I] 5577','[O_I] 6300','[S_III] 6312','[O_I] 6363','[N_II] 6548','H_alpha','[N_II] 6583','[S_II] 6716','[S_II] 6730','[Ar_III] 7135']

# The list of properties of each line in above list
SpZprops=['plate','mjd','fiberID','name','wave','z','z_e','s','s_e','a','a_e','w','w_e','c','c_e','npix','dof','chi2']
#PLATE 	Plate number
#MJD 	Modified Julian date of observation
#FIBERID 	Fiber ID (1 to 640)
#LINENAME 	Line name
#LINEWAVE 	Catalog wavelength for this line in vacuum Angstroms
#LINEZ 	Redshift
#LINEZ_ERR 	Redshift error (negative for invalid fit)
#LINESIGMA 	Gaussian width in km/sec
#LINESIGMA_ERR 	Error in gaussian width (negative for invalid fit)
#LINEAREA 	Area in gaussian fit where units are (flux-units) * Ang
#LINEAREA_ERR 	Flux error (negative for invalid fit)
#LINEEW 	Equivalent width (Angstroms)
#LINEEW_ERR 	Equivalent width error (negative for invalid fit)
#LINECONTLEVEL 	Continuum level at line center
#LINECONTLEVEL_ERR 	Error in continuum level at line center
#LINENPIX 	Number of good pixels within +/- 3 sigma of the line center
#LINEDOF 	Degrees of freedom in fit, approximated as LINENPIX minus the number of terms fit for that line, which can be fractional if one parameter if fixed betwen several lines
#LINECHI2 	χ2 for all points within +/- 3 sigma of the line center (negative if no such points)


# HELPER FUNCTIONS
TYPES={'agn':'INTEGER',
       'ObjID':'INTEGER PRIMARY KEY',
       'SpecObjID':'INTEGER',
       'flags':'INTEGER',
       'zwarning ':'INTEGER',
       'plate':'INTEGER',
       'mjd':'INTEGER',
       'fiberID':'INTEGER',
       'fit_type_old':'INTEGER',
       'fit_type_young':'INTEGER',
       'primtarget':'INTEGER'
       }
def createcolumnifnotexists(curs,name,ctype=None,table='sdss'):
    try: curs.execute('SELECT %s from %s LIMIT 1;'%(name,table))
    except Exception,e:
        if ctype: ty = ctype
        elif TYPES.has_key(name): ty=TYPES[name]
        else: ty='REAL'
        sql='ALTER TABLE %s ADD COLUMN %s %s;'%(table,name,ty)
        print sql
        curs.execute(sql)


def list2csv(list,delim=DELIM):
    return string.join(list,delim)

def printcomm():
    print 'Do not forget to commit the connection!'

def newobject(cols,line,cursor):
    id=line[cols.index('objID')]
    cursor.execute('SELECT objID from sdss WHERE objID==%s'%id)
    try:
        cursor.next()
        return False
    except: return True

def insert(cols,line,cursor):
    cursor.execute('INSERT INTO sdss (%s) VALUES (%s)'%(list2csv(cols),list2csv(line)))

def update(cols,line,cursor):
    id=line[cols.index('objID')]
    assignments=''
    for i,c in enumerate(cols):
        if (c != 'objID'):
            assignments+='%s=%s,'%(c,line[i])
    cursor.execute('UPDATE sdss SET %s WHERE objID==%s'%(assignments[:-1],id))

def setupdb(dbname=DBNAME):
    connection=sqlite.connect(dbname)
    cursor=connection.cursor()

    return connection,cursor

def createviews(cursor):
    cursor.execute('CREATE VIEW clean AS SELECT * from sdss WHERE zconf>0.95 AND ((flags & 8)==0) AND (primtarget & 64 > 0) AND z > 0.005')
    #cursor.execute('CREATE VIEW clean AS SELECT * from sdss WHERE zconf>0.95')
    cursor.execute('CREATE VIEW sb AS SELECT * from clean WHERE Ha_w > 60.0 AND Ha_s>1')
    cursor.execute('CREATE VIEW pb AS SELECT * from clean WHERE Hd_w < -6.0')
    cursor.execute('create view sball as select * from clean inner join sbfit on sbfit.ID=clean.specObjID;')
    cursor.execute('create view pball as select * from clean inner join pbfit on pbfit.ID=clean.specObjID;')
    printcomm()

def readfile(fname,cursor,table='sdss',delim=DELIM,cnames=None,ctypes=None):
    file=open(fname)
    if not cnames:
        cols=map(string.strip,file.readline().split(delim))
    if not ctypes:
        types=map(string.strip,file.readline().split(delim))
    poplist=[] # the ones we really want
    for i,c in enumerate(cols):
        if c=='kasta':
            print 'skipping column %d'%i
            poplist.append(i)
            continue
        createcolumnifnotexists(cursor,c,types[i],table=table)

    poplist=N.array(poplist)-N.arange(len(poplist))
    for p in poplist:
        cols.pop(p)

    for line in file:
        line=map(string.strip,line.split(DELIM))
        for p in poplist: line.pop(p)
        if newobject(cols,line,cursor): insert(cols,line,cursor)
        else: update(cols,line,cursor)

    file.close()

def getspecfilename(mjd,plate,fiberID,addBase=True):
    if addBase: b = SPECBASE
    else: b = ''
    return join(b,'spSpec-%05d-%04d-%03d.fit'%(mjd,plate,fiberID))

def getspecZfilename(mjd,plate):
    return join(SPECZBASE,'spZline-%04d-%05d.fits'%(plate,mjd))


def specfromid(id,cursor=False,db=DBNAME):
    if cursor==False:
        connection,cursor=setupdb(db)
        closeafter=True
    else: closeafter=False

    cursor.execute('SELECT mjd,plate,fiberID FROM sdss WHERE objID==%s'%id)
    fits=F.open(getspecfilename(cursor.next()))

    if closeafter:
        cursor.close()
        connection.close()
    return fits

def specsfromquery(query,cursor=False,db=DBNAME):
    if cursor==False:
        connection,cursor=setupdb(db)
        closeafter=True

    if closeafter:
        cursor.close()
        connection.close()
    return speclist

def fetch(cursor):
    result = map(N.array,zip(*cursor.fetchall()))
    if len(result) == 1: result = result[0]
    return result

def get(cursor,question):
    print question
    cursor.execute(question)
    return fetch(cursor)

def gettable(cursor,cols,where=None,table='sdss'):
    if where: where=' WHERE '+where
    return get(cursor,'SELECT %s FROM %s%s'%(cols,table,where))

def getsb(cursor,cols,where=None):
    return gettable(cursor,cols,where=where,table='sball')

def getpb(cursor,cols,where=None):
    return gettable(cursor,cols,where=where,table='pball')

def getspZline(cursor,wantedlines=[11,12,15,24],wantedprops=N.arange(7,15)):
    for l in wantedlines:
        for p in wantedprops:
            colname='p_%02d_%s'%(l,SpZprops[p])
            try: cursor.execute('ALTER TABLE sdss ADD COLUMN %s REAL'%colname)
            except: print 'Colums %s exists already'%colname
    cursor.execute('SELECT plate,mjd,fiberID,objID FROM sdss ORDER BY plate,mjd,fiberID')
    plate,mjd,fiberID,objID=fetch(cursor)
    for i,obj in enumerate(objID):
        data=F.open(getspecZfilename((mjd[i],plate[i])))[1].data
        data.shape=640,-1
        comm=''
        for j in wantedlines:
            dat=data[fiberID[i]-1,j]
            if (dat[0] != plate[i]) or (dat[1] != mjd[i]) or (dat[2] != fiberID[i]):
                print 'somethings terribly wrong'
                break
            for k in wantedprops:
                comm+='%s=%s,'%('p_%02d_%s'%(j,SpZprops[k]),dat[k])

        comm='UPDATE sdss SET '+comm[:-1]+' WHERE objID==%s'%objID[i]
        #print comm
        cursor.execute(comm)
    printcomm()


def dumpascispec(curs,where='z<5',table='sdss'):
    all=curs.execute('SELECT objID,mjd,plate,fiberID,Ha_h,Ha_w,Hb_h,Hb_w,z from %s WHERE %s ORDER BY objID'%(table,where))
    for objID,mjd,plate,fiberID,Ha_h,Ha_w,Hb_h,Hb_w,z in all:
        fname=getspecfilename((mjd,plate,fiberID))
        fits=F.open(fname)
        head,spec,noise=sdss.splitfits(fits)
        wave=10**(head.get('COEFF0')+(N.arange(len(spec),dtype='f')*head.get('COEFF1')))
        outf=open('%d.spec'%objID,'w')
        outf.write('%.8e\n'%z)
        outf.write('%.8e\n'%Ha_h)
        outf.write('%.8e\n'%Ha_w)
        outf.write('%.8e\n'%Hb_h)
        outf.write('%.8e\n'%Hb_w)
        for i,s in enumerate(spec):
            outf.write('%.8e %.8e\n'%(wave[i],s))
        outf.close()



#################

def main():
    print 'no main(). this is a library.'

if __name__=='__main__':
    main()
