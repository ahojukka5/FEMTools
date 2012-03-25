# -*- coding: utf-8 -*-

from __future__ import division

import teefem
import logging
import numpy as np
from common import Model
from numpy import array, zeros, matrix
from numpy import sqrt
from teefem.boundary_conditions import PressureBoundaryCondition

###############################################################################

class PlaneStressElement1D(teefem.elements.Element1D):
    '''
    Plane stress 1d elements, general class.
    '''
    
    def __init__(self, *args, **kwds):
        super(PlaneStressElement1D, self).__init__(*args, **kwds)
        self.pressure = lambda k: 0

    def assign_boundary_condition(self, bc):
        if bc.__class__ is PressureBoundaryCondition:
            self.pressure = bc.pressure

    @property
    def force_vector(self):
        Jk = self.geom.Jk
        Je = self.geom.Je
        N = self.geom.N
        Rq = zeros(self.dimension)
        qx = array([W * self.pressure(k) * N(k) * Jk(k) for (W,k) in zip(self.iweights,self.ipoints)])
        qy = array([W * self.pressure(k) * N(k) * Je(k) for (W,k) in zip(self.iweights,self.ipoints)])
        Rq[1::2] += np.sum(qy,axis=0)
        Rq[0::2] += np.sum(qx,axis=0)
        return Rq

class MEPLSE2(PlaneStressElement1D):
    ''' 
    Plane stress linear line element.
    '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*2
        super(MEPLSE2, self).__init__(*args, **kwds)
        self.ipoints = [(0.0)]
        self.iweights = [(2.0)]


class MEPLSE3(PlaneStressElement1D):
    '''
    Plane stress quadratic line element.
    '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*3
        super(MEPLSE3, self).__init__(*args, **kwds)
        self.ipoints = [(-1/sqrt(3)),(1/sqrt(3))]
        self.iweights = [(1.0),(1.0)]

###############################################################################

class PlaneStressElement2D(teefem.elements.Element2D):
    '''
    Plane stress 2d elements, general class.
    '''

    def __init__(self, *args, **kwds):
        super(PlaneStressElement2D, self).__init__(*args, **kwds)
        self.dimension = len(self.geom.nodes)*2 # (dx,dy)
        self.has_stiffness = True
        self.dNdk = self.geom.dNdk
        self.dNde = self.geom.dNde
    
    @property
    @teefem.cache
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
    
    def update(self, U):
        ''' '''
        u = zeros(self.dimension)
#        print("dim: {0}".format(u.shape))
        for i in xrange(len(self.nodes)):
            node = self.nodes[i]
            node.update_field(
                field = 'DEPL', 
                params = ('DX','DY'), 
                values = (U[node.gdof[0]], U[node.gdof[1]]))
            for j in xrange(2):
#                print i*2+j
                u[i*2+j] = U[node.gdof[j]]
        self.u = matrix(u).T

    @teefem.cache
    def S(self,*ke):
        B = self.B(*ke)
        D = self.D(*ke)
        return D*B*self.u

    @teefem.cache
    def vmis(self, *ke):
        S = self.S(*ke)
#        print S
        sxx = S[0,0]
        syy = S[1,0]
        sxy = S[2,0]
        return np.sqrt(sxx**2 + syy**2 - sxx*syy + 3*sxy**2)
        

class MECPTR3(PlaneStressElement2D):
    """ Plane stress linear triangle finite element CST
    
    EXAMPLES:     
    
    >>> n1 = Node(x=0.0, y=0.0)
    >>> n2 = Node(x=2.0, y=1.0)
    >>> n3 = Node(x=0.5, y=2.0)
    >>> tr3 = Tria3(name = 'tr3', nodes = (n1, n2, n3))
    >>> tr3.material = Elastic(E = 210e9, nu = 0.3)
    >>> el3 = MECPTR3(shape = tr3)
    >>> el3.plot(filename = 'img/MECPTR3.png')
    
    .. image:: img/MECPTR3.png

    >>> el3emat = el3.material_matrix
    >>> print el3emat
    [[  2.3e+11   6.9e+10   0.0e+00]
     [  6.9e+10   2.3e+11   0.0e+00]
     [  0.0e+00   0.0e+00   8.1e+10]]

    >>> el3bmat = el3.kinematic_matrix(1/3,1/3)
    >>> print el3bmat
    [[-0.3  0.   0.6  0.  -0.3  0. ]
     [ 0.  -0.4  0.  -0.1  0.   0.6]
     [-0.4 -0.3 -0.1  0.6  0.6 -0.3]]
     
    >>> el3kmat = el3.stiffness_matrix
    >>> print el3kmat
    [[  5.9e+10   3.2e+10  -5.7e+10  -3.0e+10  -1.6e+09  -2.5e+09]
     [  3.2e+10   8.6e+10  -2.4e+10   1.6e+09  -8.2e+09  -8.7e+10]
     [ -5.7e+10  -2.4e+10   1.3e+11  -2.1e+10  -7.7e+10   4.5e+10]
     [ -3.0e+10   1.6e+09  -2.1e+10   5.4e+10   5.1e+10  -5.6e+10]
     [ -1.6e+09  -8.2e+09  -7.7e+10   5.1e+10   7.9e+10  -4.3e+10]
     [ -2.5e+09  -8.7e+10   4.5e+10  -5.6e+10  -4.3e+10   1.4e+11]]
    """    
    
    def __init__(self, *args, **kwds):
        super(MECPTR3, self).__init__(*args, **kwds)
        self.ipoints = [(1.0/3.0, 1.0/3.0)]
        self.iweights = [(1.0/2.0)]
        
class MECPTR6(PlaneStressElement2D):
    """ Plane stress quadratic triangular finite element 
    
    EXAMPLES:     
    
    >>> n1 = Node(x=0.0, y=0.0)
    >>> n2 = Node(x=2.0, y=1.0)
    >>> n3 = Node(x=0.5, y=2.0)
    >>> n4 = Node(x=1.0, y=0.0)
    >>> n5 = Node(x=1.0, y=1.0)
    >>> n6 = Node(x=0.2, y=1.0)
    >>> tr6 = Tria6(name = 'tr6', nodes = (n1, n2, n3, n4, n5, n6))
    >>> tr6.material = Elastic(E = 210e9, nu = 0.3)
    >>> el6 = MECPTR6(shape = tr6)
    >>> el6.plot(filename = 'img/MECPTR6.png')
    
    .. image:: img/MECPTR6.png

    >>> el6emat = el6.material_matrix
    >>> print el6emat
    [[  2.3e+11   6.9e+10   0.0e+00]
     [  6.9e+10   2.3e+11   0.0e+00]
     [  0.0e+00   0.0e+00   8.1e+10]]

    >>> el6bmat = el6.kinematic_matrix(1/3,1/3)
    >>> print el6bmat
    [[-0.2  0.   0.2  0.  -0.   0.   0.1  0.   0.7  0.  -0.8  0. ]
     [ 0.  -0.2  0.  -0.   0.   0.2  0.  -0.7  0.   0.6  0.   0.1]
     [-0.2 -0.2 -0.   0.2  0.2 -0.  -0.7  0.1  0.6  0.7  0.1 -0.8]]
     
    >>> el6kmat = el6.stiffness_matrix
    >>> print el6kmat
    [[  1.9e+11   7.2e+10   2.3e+10   1.3e+10   1.9e+10   1.7e+10  -1.5e+11  -2.0e+10  -9.5e+09  -2.5e+10  -7.4e+10  -5.7e+10]
     [  7.2e+10   1.4e+11   1.1e+10   8.4e+09   1.9e+10   2.7e+10  -1.2e+10  -4.4e+10  -2.5e+10  -1.2e+09  -6.4e+10  -1.3e+11]
     [  2.3e+10   1.1e+10   1.2e+11   4.6e+09   1.5e+10  -1.0e+10  -1.5e+10  -4.4e+10  -1.5e+11   4.8e+10   1.2e+10  -9.7e+09]
     [  1.3e+10   8.4e+09   4.6e+09   4.2e+10  -1.2e+10   7.1e+09  -5.1e+10  -2.1e+10   5.6e+10  -4.1e+10  -9.7e+09   3.7e+09]
     [  1.9e+10   1.9e+10   1.5e+10  -1.2e+10   3.6e+10   9.3e+09   9.9e+08  -2.2e+10   5.3e+09   7.0e+10  -7.5e+10  -6.5e+10]
     [  1.7e+10   2.7e+10  -1.0e+10   7.1e+09   9.3e+09   8.5e+10  -2.2e+10   1.1e+10   6.2e+10  -3.6e+10  -5.7e+10  -9.5e+10]
     [ -1.5e+11  -1.2e+10  -1.5e+10  -5.1e+10   9.9e+08  -2.2e+10   2.6e+11  -5.5e+10  -6.6e+10   4.1e+10  -2.8e+10   9.9e+10]
     [ -2.0e+10  -4.4e+10  -4.4e+10  -2.1e+10  -2.2e+10   1.1e+10  -5.5e+10   3.6e+11   4.1e+10  -2.7e+11   9.9e+10  -4.5e+10]
     [ -9.5e+09  -2.5e+10  -1.5e+11   5.6e+10   5.3e+09   6.2e+10  -6.6e+10   4.1e+10   5.9e+11  -9.6e+10  -3.7e+11  -3.8e+10]
     [ -2.5e+10  -1.2e+09   4.8e+10  -4.1e+10   7.0e+10  -3.6e+10   4.1e+10  -2.7e+11  -9.6e+10   4.4e+11  -3.8e+10  -9.4e+10]
     [ -7.4e+10  -6.4e+10   1.2e+10  -9.7e+09  -7.5e+10  -5.7e+10  -2.8e+10   9.9e+10  -3.7e+11  -3.8e+10   5.4e+11   7.1e+10]
     [ -5.7e+10  -1.3e+11  -9.7e+09   3.7e+09  -6.5e+10  -9.5e+10   9.9e+10  -4.5e+10  -3.8e+10  -9.4e+10   7.1e+10   3.6e+11]]
    """
    def __init__(self, *args, **kwds):
        super(MECPTR6, self).__init__(*args, **kwds)
        self.ipoints = [(1.0/6.0, 1.0/6.0),(2.0/3.0,1.0/6.0),(1.0/6.0,2.0/3.0)]
        self.iweights = [1.0/6.0, 1.0/6.0, 1.0/6.0]

class MECPQU4(PlaneStressElement2D):
    """ Plane stress linear quadrangle finite element
    
    EXAMPLES:     
    
    >>> n1 = Node(x=1.0, y=0.0)
    >>> n2 = Node(x=2.5, y=0.5)
    >>> n3 = Node(x=2.0, y=2.0)
    >>> n4 = Node(x=-0.5, y=1.5)
    >>> qu4 = Quad4(name = 'MA1', nodes = (n1, n2, n3, n4))
    >>> qu4.material = Elastic(E = 100e9, nu = 0.5)
    >>> ma1 = MECPQU4(shape = qu4)
    >>> ma1.plot(filename = 'img/MECPQU4.png')
    
    .. image:: img/MECPQU4.png

    >>> ma1emat = ma1.material_matrix
    >>> print ma1emat
    [[  1.3e+11   6.7e+10   0.0e+00]
     [  6.7e+10   1.3e+11   0.0e+00]
     [  0.0e+00   0.0e+00   3.3e+10]]

    >>> ma1bmat = ma1.kinematic_matrix(0,0)
    >>> print ma1bmat
    [[-0.1  0.   0.3  0.   0.1  0.  -0.3  0. ]
     [ 0.  -0.4  0.  -0.1  0.   0.4  0.   0.1]
     [-0.4 -0.1 -0.1  0.3  0.4  0.1  0.1 -0.3]]
     
    >>> ma1kmat = ma1.stiffness_matrix
    >>> print ma1kmat
    [[  4.7e+10   2.3e+10  -3.0e+10  -1.1e+10  -1.9e+10  -2.0e+10   2.0e+09   8.6e+09]
     [  2.3e+10   1.1e+11  -2.8e+10  -2.7e+09  -2.0e+10  -7.0e+10   2.5e+10  -3.9e+10]
     [ -3.0e+10  -2.8e+10   6.1e+10  -1.2e+10  -1.5e+09   2.5e+10  -2.9e+10   1.5e+10]
     [ -1.1e+10  -2.7e+09  -1.2e+10   4.9e+10   8.3e+09  -4.4e+10   1.5e+10  -2.5e+09]
     [ -1.9e+10  -2.0e+10  -1.5e+09   8.3e+09   4.0e+10   2.2e+10  -1.9e+10  -1.0e+10]
     [ -2.0e+10  -7.0e+10   2.5e+10  -4.4e+10   2.2e+10   1.0e+11  -2.7e+10   1.3e+10]
     [  2.0e+09   2.5e+10  -2.9e+10   1.5e+10  -1.9e+10  -2.7e+10   4.7e+10  -1.4e+10]
     [  8.6e+09  -3.9e+10   1.5e+10  -2.5e+09  -1.0e+10   1.3e+10  -1.4e+10   2.8e+10]]
    
    """
    def __init__(self, *args, **kwds):
        super(MECPQU4, self).__init__(*args, **kwds)
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf s. 10
        a = 1/sqrt(3)
        self.ipoints = [(-a, -a),(a,-a),(a,a),(-a,a)]
        self.iweights = [(1.0),(1.0),(1.0),(1.0)]

class MECPQU8(PlaneStressElement2D):
    ''' Plane stress quadratic 8-node Serendip finite element '''
    def __init__(self, *args, **kwds):
        super(MECPQU8, self).__init__(*args, **kwds)
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf s. 10
        a = sqrt(3/5)
        self.ipoints = [(-a, -a),(a,-a),(a,a),(-a,a),(0,-a),(a,0),(0,a),(-a,0),(0,0)] # Integroimispisteet
        self.iweights = [(25/81),(25/81),(25/81),(25/81),(40/81),(40/81),(40/81),(40/81),(64/81)] # Painokertoimet

class MECPQU9(PlaneStressElement2D):
    ''' Plane stress quadratic 9-node finite element '''
    def __init__(self, *args, **kwds):
        super(MECPQU9, self).__init__(*args, **kwds)
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf s. 10
        a = sqrt(3/5)
        self.ipoints = [(-a, -a),(a,-a),(a,a),(-a,a),(0,-a),(a,0),(0,a),(-a,0),(0,0)] # Integroimispisteet
        self.iweights = [(25/81),(25/81),(25/81),(25/81),(40/81),(40/81),(40/81),(40/81),(64/81)] # Painokertoimet

###############################################################################

class C_PLAN(Model):
    
    ''' Plane stress modelisation '''
    
    def __init__(self, *args, **kwds):
        super(C_PLAN, self).__init__(*args, **kwds)
        self.mapping = {
            'Seg2': MEPLSE2,
            'Seg3': MEPLSE3,
            'Tria3': MECPTR3,
            'Tria6': MECPTR6,
            'Quad4': MECPQU4,
            'Quad8': MECPQU8,
            'Quad9': MECPQU9,
        }
        self.nodedofs = ('dx','dy')
        self.nodeloads = ('fx','fy')
        self.nodedim = len(self.nodedofs)

        self.init()
    
cplan = C_PLAN

logging.debug("Module {0} loaded.".format(__file__))