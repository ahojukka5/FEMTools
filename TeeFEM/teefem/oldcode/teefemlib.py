# -*- coding: utf-8 -*-
"""
Created on Wed Feb 01 19:43:07 2012

@author: Jukka Aho <jukka.aho@kapsi.fi>

TeekkariFEM

Ratkaisee tasojännitystilan probleemia.

http://assets.en.oreilly.com/1/event/27/Best%20practices%20for%20_scripting_%20with%20Python%203%20Paper.pdf

"""

from __future__ import division
# import tables
from numpy import sqrt, zeros, array, matrix, sum, linspace, concatenate
from numpy.linalg import det, inv
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
import sys
import os
import locale
import numpy as np
import unittest
import re
import matplotlib.pylab as plt

np.set_printoptions(precision=1,linewidth=160)
np.set_printoptions(edgeitems=100)
enc = locale.getpreferredencoding()



import materials

#####################################
### Helpers.py - HELPER FUNCTIONS ###
#####################################

def dist(n1,n2):
    ''' Returns distance between two nodes '''
    return sqrt((n2.x-n1.x)**2 + (n2.y-n1.y)**2 + (n2.z-n1.z)**2)

def arrtostr(arr):
    ''' Converts list of integers to string '''
    return ''.join([chr(i) for i in arr]).strip()

class IJV(object):
    """ Creates IJV-arrays and converts to COO sparse matrix 
    
    INPUT:
        - ``dim`` - lenght of I,J,V lists
    
    EXAMPLES:
        
    >>> k = 210e9*10e-4/3
    >>> ijv = IJV(4)
    >>> ijv.add(0,0,k)
    >>> ijv.add(1,1,k)
    >>> ijv.add(0,1,-k)
    >>> ijv.add(1,0,-k)
    >>> print ijv.tocoo().todense()
    [[ 70000000. -70000000.]
     [-70000000.  70000000.]]
    
    """
    def __init__(self,dim):
        self.idx = 0
        self.dim = dim
        self.i = zeros(dim,dtype=int)
        self.j = zeros(dim,dtype=int)
        self.v = zeros(dim)
        
    def add(self,i,j,v):
        self.i[self.idx] = i
        self.j[self.idx] = j
        self.v[self.idx] = v
        self.idx += 1
        
    def tocoo(self):
        return coo_matrix((self.v,(self.i,self.j)))


##########################
### 1D finite elements ###
##########################


##########################
### Models.py - MODELS ###
##########################

class Model(object):
    ''' General model function '''
    def __init__(self, *args, **kwds):
        self.mesh = kwds['mesh']
        self.elements = {}
        # Isoparametric mappign -> same nodes
        self.nodes = self.mesh.nodes
        self.group_no = self.mesh.group_no
        logging.debug("Model created")

    def assign_bc(self,*args,**kwds):
        ''' Assigns boundary conditions to elements / nodes '''
        if kwds.has_key('elem_group'):
            for element in self.group_ma[kwds['elem_group']]:
                element.bc = kwds
        if kwds.has_key('node_group'):
            try:
                for node in self.group_no[kwds['node_group']]:
                    node.bc = kwds
            except KeyError:
                logging.debug("Error: node group {0} does not exist!".format(kwds['node_group']))


#############################
### Tests.py - UNIT TESTS ###
#############################

class C_PLAN_tests(unittest.TestCase):
    
    """ Unit test for plane stress modelisation.
        Plate bend, results are compared to Code Aster results. 
    """
    
    def runTest(self):
        
        data = (
        {'testname': 'Linear triangle elements MECPTR3',
         'meshfile': 'mesh1.mail',
         'refresult': -1.00963370103376E-02
         },
        {'testname': 'Quadratic triangle elements MECPTR6',
         'meshfile': 'mesh2.mail',
         'refresult': -1.06681040183457E-02
         },
        {'testname': 'Linear quadrangle elements MECPQU4',
         'meshfile': 'mesh3.mail',
         'refresult': -1.06819027007511E-02
         },
        {'testname': 'Quadratic quadrangle elements MECPQU8',
         'meshfile': 'mesh4.mail',
         'refresult': -1.07998813565835E-02
         },
        )
        
        for d in data:
            print(d['testname'])
            E = 210e9
            nu = 0.3
            t = 10e-3
            F = 100e3
            q = 100e3
            mesh = Mesh(filename = d['meshfile'])
            model = C_PLAN(mesh = mesh)
            mat = Elastic(E = E, nu = nu)
            mesh.assign_material(mat)
            # Huom. t = 1 !!!!
            model.assign_bc(elem_group = 'GA2', pressure = q/t)
            model.assign_bc(node_group = 'P1', fy = -F/t)
            model.assign_bc(node_group = 'GA1', dx = 0, dy = 0)
            #model.static_solve()
            K = model.assembly_stiffness_matrix().tocsr()
            R = model.assembly_force_vector()
            U = spsolve(K,R)
            for node in model.nodes.itervalues():
                node.update_field(
                    field = 'DEPL', 
                    params = ('DX','DY'), 
                    values = (U[node.gdof[0]], U[node.gdof[1]]),
                    )
            #print("Displacement DY on P1")
            P1 = iter(model.group_no['P1']).next()
            dy = P1.fields['DEPL']['DY']
            CA = d['refresult'] # Code Asterin ratkaisu
            print("Solution (Python): %0.14E"%(dy))
            print("Solution (Aster) : %0.14E"%(CA))
            diff = abs(dy-CA)
            print("Diff : %0.14E"%(diff))
            eps = 1e-10
            self.assertTrue(diff<eps)

#################
### MAIN INIT ###
#################

def main():
    ''' Main program '''
    print("Running tests")
    import doctest
    doctest.testmod()
    unittest.main()
    print("Tests done.")

if __name__ == '__main__':
    main()
