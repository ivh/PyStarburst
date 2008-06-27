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

# GLOBAL VARIABLES
DBNAME=join('/','data','sdss','sb.db')
DELIM=','
SPECBASE=join('/','data','sdss','spectra')

# HELPER FUNCTIONS
def list2csv(list):
    return string.join(list,DELIM)

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
    if not exists(dbname): newdb=True
    else: newdb=False

    connection=sqlite.connect(dbname)
    cursor=connection.cursor()

    if newdb:
        cursor.execute('CREATE TABLE sdss (objID INTEGER PRIMARY KEY)')
        connection.commit()

    return connection,cursor

def readfiles(fnames,cursor):
    for fname in fnames:
        print "Working on file: %s"%fname
        file=open(fname)
        cols=map(string.strip,file.readline().split(DELIM))
        types=map(string.strip,file.readline().split(DELIM))
        print 'found columns: %s'%cols
        for i,c in enumerate(cols):
            comm='ALTER TABLE sdss ADD COLUMN %s %s'%(c,types[i])
            try: cursor.execute(comm)
            except: print 'not adding pre-existing column: %s'%c
        
        for line in file:
            line=map(string.strip,line.split(DELIM))
            if newobject(cols,line,cursor): insert(cols,line,cursor)
            else: update(cols,line,cursor)
        
        file.close()

def getspecfilename((mjd,plate,fiberID)):
    return join(SPECBASE,'spSpec-%05d-%04d-%03d.fit'%(mjd,plate,fiberID))
    

def specfromid(id,cursor=False,db=DBNAME):
    if cursor==False:
        connection,cursor=setupdb(db)
        closeafter=True

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


#################

def usage_example():
    connection,cursor=setupdb()
    cursor.execute('SELECT Ha_w,Hb_w FROM sdss WHERE zConf>0.95')
    b=N.array(map(list,cursor.fetchall()))
    print b.shape
    P.plot(b[:,0],b[:,1],'ro')
    P.show()


def main():
  
    if sys.argv[1].endswith('.db'):
        connection,cursor=setupdb(sys.argv[1])
        readfiles(sys.argv[2:],cursor)
    else:
        connection,cursor=setupdb()
        readfiles(sys.argv[1:],cursor)

    connection.commit()
    cursor.close()
    connection.close()

    #usage_example()


if __name__=='__main__':
    main()



# query stuff
#cursor.execute('SELECT email,acctype FROM participants')
#print cursor.fetchall()
#names=cursor.fetchmany()
#names=cursor.fetchone()
#names=cursor.next()
