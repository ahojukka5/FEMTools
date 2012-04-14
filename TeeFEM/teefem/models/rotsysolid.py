# -*- coding: utf-8 -*-

# Rotsysolid

from __future__ import division

import teefem
import logging
import numpy as np
from common import Model
from numpy import zeros, matrix

class RotSolid(teefem.elements.Element2D):
    '''
    Rotsy solid
    '''

    def __init__(self, *args, **kwds):
        super(RotSolid, self).__init__(*args, **kwds)
        self.dimension = len(self.geom.nodes)*2 # (dx,dy)
        self.degrees_of_freedom = ('DX','DY')
        self.has_stiffness = True
        self.dNdk = self.geom.dNdk
        self.dNde = self.geom.dNde
        self.ipoints = [(1.0/3.0, 1.0/3.0)]
        self.iweights = [(1.0/2.0)]    
    
    @property
    def detJ(self):
        return self.geom.detJ

    def assign_material(self, mat):
        self.material = mat

    @teefem.cache    
    def B(self,*ke):
        """
        Returns plane kinematic matrix B 
        """
        dNdk = self.dNdk(*ke)
        dNde = self.dNde(*ke)
        dNxy = self.invJ(*ke) * matrix([dNdk,dNde])
        dNdx = dNxy[0]
        dNdy = dNxy[1]
        B = matrix(zeros((3,self.dimension)))
        B[0,0::2] = dNdx
        B[1,1::2] = dNdy
        B[2,0::2] = dNdy
        B[2,1::2] = dNdx
        return B

    @teefem.cache
    def D(self,*ke):
        """
        Material matrix D
        """
        EE = self.material.E
        nu = self.material.nu
        return EE/(1-nu**2)*matrix([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])

    @property
#    @teefem.cache
    def stiffness_matrix(self):
        ''' Stiffness matrix K '''
        K = matrix(zeros((self.dimension,self.dimension)))
        for (W,ke) in zip(self.iweights, self.ipoints):
            detJ = self.detJ(*ke)
            B = self.B(*ke)
            D = self.D(*ke)
            K += W*B.T*D*B*detJ
        return K

class ROTSOL(Model):
    
    ''' Plane stress modelisation '''
    
    def __init__(self, *args, **kwds):
        super(ROTSOL, self).__init__(*args, **kwds)
        self.mapping = {
            'Tria3': RotSolid,
#            'Tria6': MECPTR6,
        }
        self.nodedofs = ('dx','dy')
        self.nodeloads = ('fx','fy')
        self.nodedim = len(self.nodedofs)

        self.init()
    
logging.debug("Module {0} loaded.".format(__file__))

###############################################################################
###############################################################################

def test1():
    ''' Palauttaa elementin jäykkyysmatriisin '''
    n1 = teefem.geom.Node(x=0, y=0, z=0)
    n2 = teefem.geom.Node(x=1, y=0, z=0)
    n3 = teefem.geom.Node(x=0, y=1, z=0)
    shape1 = teefem.geom.Tria3(nodes=(n1,n2,n3))
    print shape1
    print dir(shape1)
    print shape1.nodes
    print shape1.detJ(0.1,0.1)
    print shape1.dNde(0.1,0.1)
    element = RotSolid(geom = shape1)
    print element
    print dir(element)
    print element.status
    print element.B(0.1,0.1)
    print element.iweights
    print element.ipoints
    
    # Tarvitaan materiaalimäärittely jäykkyysmatriisia varten,
    # kokeile commentata -> virheilmoitus
    element.material = teefem.materials.elastic(E=210e9, nu=0.3)

    print element.stiffness_matrix
    
    import matplotlib.pylab as plt
    shape1.plot()
    plt.show()

if __name__ == '__main__':
    test1()