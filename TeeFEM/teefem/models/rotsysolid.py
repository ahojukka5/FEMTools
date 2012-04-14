# -*- coding: utf-8 -*-
"""

Muutokset

14.4.2012
- Lisätty RotSolid __init__: has_stiffness_matrix , has_load_vector
- Lisätty yksinkertainen piirtorutiini
Jukka

"""

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
        self.has_stiffness_matrix = True
        self.has_load_vector = False
        self.dNdk = self.geom.dNdk
        self.dNde = self.geom.dNde
        self.N = self.geom.N
        self.x = self.geom.x
        self.ipoints = [[1.0/6.0, 1.0/6.0],[2.0/3.0, 1.0/6.0], [1.0/6.0, 2.0/3.0]]
        self.iweights = [[1.0/3.0],[1.0/3.0],[1.0/3.0]]    
    
    @property
    def detJ(self):
        return self.geom.detJ

    def assign_material(self, mat):
        self.material = mat
    
    def B(self,*ke):
        """
        Returns plane kinematic matrix B 
        """
        dNdk = self.dNdk(*ke)
        dNde = self.dNde(*ke)
        N = self.N(*ke)
        r = self.x(*ke)
        dNxy = self.invJ(*ke) * matrix([dNdk,dNde])
        dNdx = dNxy[0]
        dNdy = dNxy[1]
        B = matrix(zeros((self.D(*ke).shape[0],self.dimension)))
        B[0,0::2] = dNdx
        B[1,0::2] = N/r
        B[2,1::2] = dNdy
        B[3,0::2] = dNdy
        B[3,1::2] = dNdx
        return B

    #@teefem.cache
    def D(self,*ke):
        """
        Material matrix D
        """
        EE = self.material.E
        nu = self.material.nu
        return EE/((1+nu)*(1-2*nu))*matrix([[1-nu,nu,nu,0],[nu,1-nu, nu,0],[nu, nu, 1-nu, 0],[0,0,0,(1-2*nu)/2]])

    @property
#    @teefem.cache
    def stiffness_matrix(self):
        ''' Stiffness matrix K '''
        K = matrix(zeros((self.dimension,self.dimension)))
        for (W,ke) in zip(self.iweights, self.ipoints):
            detJ = self.detJ(*ke)
            r=self.x(*ke)
            B = self.B(*ke)
            D = self.D(*ke)
            K += np.pi*W[0]*B.T*D*B*r*detJ
        return K

class ROTSOL(Model):
    
    ''' Plane stress modelisation '''
    
    def __init__(self, *args, **kwds):
        super(ROTSOL, self).__init__(*args, **kwds)
        self.mapping = {
            'Tria3': RotSolid,
        }
        self.nodedofs = ('dx','dy')
        self.nodeloads = ('fx','fy')
        self.nodedim = len(self.nodedofs)

        self.init()

    def plot(self, **kwds):
        ''' Plottaus. '''
        import matplotlib.pylab as plt
        fig = plt.figure()
        ax = fig.gca()
        ax.grid()
        sf = kwds.get('scalefactor',1)
        for e in self.stiff_elements:
            x1 = e.nodes[0].x + e.nodes[0].fields['DEPL']['DX']*sf
            y1 = e.nodes[0].y + e.nodes[0].fields['DEPL']['DY']*sf
            x2 = e.nodes[1].x + e.nodes[1].fields['DEPL']['DX']*sf
            y2 = e.nodes[1].y + e.nodes[1].fields['DEPL']['DY']*sf
            x3 = e.nodes[2].x + e.nodes[2].fields['DEPL']['DX']*sf
            y3 = e.nodes[2].y + e.nodes[2].fields['DEPL']['DY']*sf
            ax.plot([x1,x2,x3,x1],[y1,y2,y3,y1],'--ko')
        return fig,ax

#logging.debug("Module {0} loaded.".format(__file__))

###############################################################################
###############################################################################

def test1():
    ''' Palauttaa elementin jäykkyysmatriisin '''
    n1 = teefem.geom.Node(x=0, y=0, z=0)
    n2 = teefem.geom.Node(x=1, y=0, z=0)
    n3 = teefem.geom.Node(x=0, y=1, z=0)
    shape1 = teefem.geom.Tria3(nodes=(n1,n2,n3))
    #print shape1
    #print dir(shape1)
    #print shape1.nodes
    #print shape1.detJ(0.1,0.1)
    #print shape1.dNde(0.1,0.1)
    #print shape1.x(0.1,0.1)
    element = RotSolid(geom = shape1)
    element.material = teefem.materials.elastic(E=210e9, nu=0.3)
    #print element
    #print dir(element)
    #print element.status
    #print element.iweights
    #print element.ipoints
    #print element.has_stiffness
    print element.stiffness_matrix
    
    import matplotlib.pylab as plt
    shape1.plot()
    plt.show()

if __name__ == '__main__':
    test1()
