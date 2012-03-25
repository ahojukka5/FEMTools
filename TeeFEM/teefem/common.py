# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 02:04:59 2012

@author: Jukka

"""

from __future__ import division

from numpy import sqrt, zeros, array,arctan2
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve as solve
import os
import logging
from logging import handlers

# Logfile
LOGFILE = os.path.expanduser('~/.teefem.log')

# create logger
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

# create formatter
#log_format = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
log_format = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# create rotating file handler
fhandler = handlers.RotatingFileHandler(LOGFILE, maxBytes=500, backupCount=2)
fhandler.setLevel(logging.DEBUG)
fhandler.setFormatter(log_format)

# crate console handler
chandler = logging.StreamHandler()
chandler.setLevel(logging.DEBUG)
chandler.setFormatter(log_format)

# add handlers to logger
#log.addHandler(fhandler)
log.addHandler(chandler)

def dist(n1,n2):
    ''' Returns distance between two nodes '''
    return sqrt((n2.x-n1.x)**2 + (n2.y-n1.y)**2 + (n2.z-n1.z)**2)

def anglxy(n1,n2):
    return arctan2(n2.y-n1.y,n2.x-n1.x)

def arrtostr(arr):
    ''' Converts list of integers to string '''
    return ''.join([chr(i) for i in arr]).strip()

class IJV(object):
    """ Creates IJV-arrays and converts to COO sparse matrix 
    
    INPUT:
        - ``dim`` - lenght of I,J,V lists
    
    EXAMPLES:
        
    >>> k = 210e9*10e-4/3
    >>> ijv = IJV()
    >>> ijv.add(0,0,k)
    >>> ijv.add(1,1,k)
    >>> ijv.add(0,1,-k)
    >>> ijv.add(1,0,-k)
    >>> print ijv.tocoo().todense()
    [[ 70000000. -70000000.]
     [-70000000.  70000000.]]
    
    """
    def __init__(self):
        log.info('Creating IJV matrix')
        self.idx = 0
        self.i = {}#zeros(dim,dtype=int)
        self.j = {}#zeros(dim,dtype=int)
        self.v = {}#zeros(dim)
        
    def add(self,i,j,v):
        self.i[self.idx] = i
        self.j[self.idx] = j
        self.v[self.idx] = v
        self.idx += 1
        
    def tocoo(self):
        log.info('Converting to COO matrix')
        v = array([self.v[x] for x in xrange(self.idx)])
        i = array([self.i[x] for x in xrange(self.idx)])
        j = array([self.j[x] for x in xrange(self.idx)])
        return coo_matrix((v,(i,j)))

def test1():
    ijv = IJV()
    ijv.add(0,0,1)
    ijv.add(1,1,1)
    ijv.add(0,1,-1)
    ijv.add(1,0,-1)
    print ijv.tocoo().todense()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
#    test1()

log.info("Module {0} loaded.".format(__file__))