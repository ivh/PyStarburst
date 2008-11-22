"""
Doc-String for the __init__.py file.
look into the doc string for the PyCigale submodule 
for more information.
"""

__author__ = 'T. Marquart'
__version__= 0.1

__PAR__= {}

import sdss as S
import db as DB
import plotting as P
import sqlcl as sqlcl
import pickle

def dump(data,filename):
  file=open(filename,'w')
  pickle.dump(data,file)
  file.close()

def load(filename):
  file=open(filename,'r')
  data=pickle.load(file)
  file.close()
  return data
