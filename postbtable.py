#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Help: none, you're on your own.

"""

import pylab as P
import numpy as N
import pyfits as F
from db import *
from numpy.ma import masked_where
from numpy import pi
from sdss import *

CRIT = 'z between 0.02 and 0.03 and b_param > 3 and Massfrac > 0.05 and agn=0 and dec < 20'

def selection2html(curs,outfile='pdbs.html',where=CRIT,table='pball',order='z asc'):
    urlbase='http://cas.sdss.org/dr7/en/tools/explore/obj.asp?id='
    wanted='objID, SpecObjID, ra, dec, z, Massfrac, b_param, E_BV, Age, ML_young, ML_tot, mass, mtot, dynMassDisk, Ha_w, Ha_s, Ha_h,Hb_h, Hd_w, m_g, m_r, f_g, f_r, Mr'
    wants=wanted.split(',')
    f=open(outfile,'w')
    f.write('<html><p>%s</p>\n'%CRIT)
    f.write('<table border=1><tr><th>#</th>\n')
    for w in wants: f.write('<th>%s</th>'%w)
    f.write('<th>M_g</th><th>M_r</th>')
    f.write('</tr>\n')
    curs.execute('SELECT DISTINCT %s FROM %s WHERE %s ORDER BY %s'%(wanted,table,where,order))
    data=curs.fetchall()
    prevID,counter=0,1
    for entry in data:
        exec(wanted+'=entry')
        f.write('<tr><td>%s</td>'%(counter,))
        if objID!=prevID: f.write('<td><a href="%s">%s</a></td>'%(urlbase+str(objID),str(objID)))
        else: f.write('<td></td>')
        f.write('<td>%s</td><td>%.5f</td><td>%.5f</td>'%(SpecObjID,ra,dec))
        f.write('<td>%.4f</td>'%z)
        f.write('<td>%.3f</td><td>%.1f</td>'%(Massfrac,b_param))
        f.write('<td>%.2f</td><td>%.1e</td><td>%.2f</td>'%( E_BV, Age, ML_young))
        f.write('<td>%.1f</td><td>%.1e</td><td>%.1e</td><td>%.1e</td>'%( ML_tot, mass, mtot, dynMassDisk))
        f.write('<td>%.1f</td><td>%.1f</td><td>%.1f</td>'%(Ha_w,Ha_s,Ha_h))
        f.write('<td>%.1f</td><td>%.1f</td>'%(Hb_h,Hd_w))
        f.write('<td>%.1f</td><td>%.1f</td>'%(m_g,m_r))
        f.write('<td>%.1f</td><td>%.1f</td><td>%.1f</td>'%(f_g,f_r,Mr))
        f.write('<td>%.1f</td>'%(absmag(m_g,z)))
        f.write('<td>%.1f</td>'%(absmag(m_r,z)))

        f.write('</tr>\n')
        prevID=objID
        counter+=1
    f.write('</table></html>')
    f.close()
