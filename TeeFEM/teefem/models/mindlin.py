# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 12:04:59 2012

@author: Jukka

Mindlin elementti, perusversio.

"""

from __future__ import division

import teefem
import teefem.elements
from teefem.boundary_conditions import PressureBoundaryCondition
from teefem.plate_functions import PlateCharacteristic
from common import Model
from teefem import log, cache

#from scipy.sparse.linalg import spsolve
from numpy import array, zeros, matrix
from numpy import sqrt
import numpy as np

class MINTR3(teefem.elements.Element2D):
    ''' MINTR3 TRIA3 Reissner-Mindlin finite element '''
    def __init__(self, *args, **kwds):
        super(MINTR3, self).__init__(*args, **kwds)
        self.dimension = 3*3
#        self.thickness = lambda ke: 0.5
#        self.has_stiffness = True
        self.kappa = 5/6
        self.degrees_of_freedom = ('DZ','DRX','DRY')

        # Kolmen pisteen integroimiskaava
        self.ipoints = [(0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]
        self.iweights = [(1.0/6.0), (1.0/6.0), (1.0/6.0)]

        # Kolmen pisteen integroimiskaava (Hammer)
#        k = 1.0/6.0
#        self.ipoints = [(k, k), (2*k, k), (k, 2*k)]
#        self.iweights = [(k), (k), (k)]

        # Interpoloidaan tuntematonta kenttää samoilla C0-jatkuvilla muoto-
        # funktioilla kuin elementin geometriaa, tässä tapauksessa siis TRIA3
        # muofofunktioilla jotka tulevat suoraan geometrialuokasta
        self.dNdk = self.geom.dNdk
        self.dNde = self.geom.dNde
        self.dNxy = lambda *ke: self.geom.invJ(*ke) * matrix([self.dNdk(*ke),self.dNde(*ke)])
        self.Nxy = self.geom.N
        self.detJ = self.geom.detJ

    def Db(self,*ke):
        E = self.material.E
        nu = self.material.nu
        return E*self.thickness(*ke)**3/(12*(1-nu**2)) * matrix([[1,nu,0],[nu,1,0],[0,0,0.5*(1-nu)]])

    def Ds(self,*ke):
        G = self.material.G
        return self.kappa * G * self.thickness(*ke) * matrix([[1,0],[0,1]])

    def assign_material(self, mat):
        self.material = mat

    def assign_boundary_condition(self, bc):
        if bc.__class__ is PressureBoundaryCondition:
            self.pressure = bc.pressure

    def assign_char(self, char):
        if char.__class__ is PlateCharacteristic:
            self.thickness = char.thickness

    def Bb(self,*ke):
        ''' Kinemaattinen matriisi Bb '''
        dNxy = self.dNxy(*ke)
        dNdx = dNxy[0]
        dNdy = dNxy[1]
        Bb = matrix(zeros((3,self.dimension)))
        Bb[0,2::3] = dNdx
        Bb[1,1::3] = dNdy
        Bb[2,1::3] = dNdx
        Bb[2,2::3] = dNdy
        return Bb

    def Bs(self,*ke):
        ''' Kinemaattinen matriisi Bs '''
        dNxy = self.dNxy(*ke)
        dNdx = dNxy[0]
        dNdy = dNxy[1]
        Nxy = self.Nxy(*ke)
        Bs = matrix(zeros((2,self.dimension)))
        Bs[0,0::3] = dNdx
        Bs[1,0::3] = dNdy
        Bs[0,2::3] = -Nxy
        Bs[1,1::3] = -Nxy
        return Bs

    @property
    def stiffness_matrix(self):
        ''' Jäykkyysmatriisi '''
        Kb = matrix(zeros((self.dimension,self.dimension)))
        Ks = matrix(zeros(Kb.shape))
        for (W,ke) in zip(self.iweights, self.ipoints):
            detJ = self.detJ(*ke)
            Bb = self.Bb(*ke)
            Db = self.Db(*ke)
            Kb += W*Bb.T*Db*Bb*detJ
            Bs = self.Bs(*ke)
            Ds = self.Ds(*ke)
            Ks += W*Bs.T*Ds*Bs*detJ
        return Kb + Ks

    @property
    def force_vector(self):
        ''' Palauttaa kuormitusvektorin '''
        b = array([W*self.pressure(*ke)*self.Nxy(*ke)*self.detJ(*ke) for (W,ke) in zip(self.iweights, self.ipoints)])
        F = zeros(9)
        F[0::3] = np.sum(b,axis=0)
        return matrix(F).T
    
    def update(self, U):
        u = zeros(self.dimension)
        nodedim = len(self.degrees_of_freedom)
        for i in xrange(len(self.nodes)):
            node = self.nodes[i]
            node.update_field(
                field = 'DEPL', 
                params = self.degrees_of_freedom, 
                values = [U[node.gdof[j]] for j in xrange(nodedim)]
                )
            for j in xrange(nodedim):
                u[i*nodedim+j] = U[node.gdof[j]]
        self.u = matrix(u).T

class MINTR3R(MINTR3):
    ''' MINTR3 TRIA3 Reissner-Mindlin finite element, reduced integration '''
    def __init__(self, *args, **kwds):
        super(MINTR3R, self).__init__(*args, **kwds)
        self.ipoints = [(1.0/3.0, 1.0/3.0)]
        self.iweights = [(1.0/2.0)]

class MindlinModel(Model):
    ''' Mindlin model '''

    def __init__(self, *args, **kwds):
        super(MindlinModel, self).__init__(*args, **kwds)
        self.nodedim = 3
        self.nodedofs = ('dz','drx','dry')
        self.nodeloads = ('fz','mx','my')

    def plot(self, **kwds):
        ''' Plottaa siirtymäkentän matplotlibillä '''
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.grid()
        for e in self.elements:
            x1 = e.nodes[0].x
            y1 = e.nodes[0].y
            z1 = e.nodes[0].z + e.nodes[0].fields['DEPL']['DZ']
            x2 = e.nodes[1].x
            y2 = e.nodes[1].y
            z2 = e.nodes[1].z + e.nodes[1].fields['DEPL']['DZ']
            x3 = e.nodes[2].x
            y3 = e.nodes[2].y
            z3 = e.nodes[2].z + e.nodes[2].fields['DEPL']['DZ']
            ax.plot([x1,x2,x3,x1],[y1,y2,y3,y1],[z1,z2,z3,z1],'--ko')
        return fig,ax

class MIN(MindlinModel):
    ''' Reissner-Mindlin model '''
    def __init__(self, *args, **kwds):
        super(MIN, self).__init__(*args, **kwds)
        self.mapping = {
            'Tria3': MINTR3,
        }
        self.init()

class MINR(MindlinModel):
    ''' Reissner-Mindlin model, reduced integration '''
    def __init__(self, *args, **kwds):
        super(MINR, self).__init__(*args, **kwds)
        self.mapping = {
            'Tria3': MINTR3R,
        }
        self.init()

mindlin = MIN
mindlin_reduced = MINR

###############################################################################
###############################################################################
###############################################################################

def test():
    
    ''' MINTR3SI-elementin matriisien tarkistus '''
    
    n1 = teefem.geom.Node(x=0, y=0, z=0.0)
    n2 = teefem.geom.Node(x=1, y=0, z=0.0)
    n3 = teefem.geom.Node(x=0, y=1, z=0.0)

    g1 = teefem.geom.Tria3(nodes = (n1,n2,n3))

    e1 = MINTR3(geom = g1)
    e1.material = teefem.materials.Elastic(E = 100.0e9, nu = 0.3)
    e1.pressure = lambda *ke: 100e3
    e1.thickness = lambda *ke: 10e-3

    print e1.status

#    import matplotlib.pylab as plt
#    g1.plot3d()

    # Aiheutaan jokin satunnainen siirtymätila
#    n1.update_field(field = 'DEPL', params = ('DZ','DRX','DRY'), values = (0.0,0.0,1.0))
#    n2.update_field(field = 'DEPL', params = ('DZ','DRX','DRY'), values = (0.0,1.0,0.0))
#    n3.update_field(field = 'DEPL', params = ('DZ','DRX','DRY'), values = (1.0,0.0,0.0))

    print("Bb(0.5,0.5):")
    print e1.Bb(0.5,0.5)
    print("Bs(0.5,0.5):")
    print e1.Bs(0.5,0.5)
    print("k:")
    k1 = e1.stiffness_matrix
    print k1
    print("sum(k.T - k): {0}".format(np.sum(k1.T - k1)))
    print("det(k)      : {0}".format(np.linalg.det(k1)))
    fq1 = e1.force_vector
    print("fq:")
    print fq1

def ex1():
    
    ''' Mindlin esimerkki 1, pyörähdyssymmetrinen laatta jakautuneella kuormalla '''
    
    t = 25e-3
    q = 100e3
    nu = 0.3
    E = 100.0e9
    a = 1
    D = E*t**3/(12*(1-nu**2))
    phi = 16/5*(t/a)**2/(1-nu)
    print phi
    w0 = q*a**4/(64*D)*(1+phi)


    mesh = teefem.mesh.unitcircle(R=1)
    
    mat = teefem.materials.elastic(E = E, nu = nu)
    mdl = MINR(mesh = mesh)
#    mdl = MIN(mesh = mesh)
    
    OM_el = mdl.elset['OM1']
    GA1_no = mdl.nset['GA1']

    teefem.assign_material(elements = OM_el, material = mat)
    
    bc1 = teefem.pressure_bc(pressure = lambda k,e: -q)
    bc2 = teefem.dirichlet_bc(encastre = True)
    
    teefem.assign_bc(elements = OM_el, bc = bc1)
    teefem.assign_bc(nodes = GA1_no, bc = bc2)

    carel = teefem.plate_functions.PlateCharacteristic(thickness = lambda k,e: t)
    teefem.assign_char(elements = OM_el, char = carel)

    mdl.static_solve()
    
    dz = [node.fields['DEPL']['DZ'] for node in mdl.nodes]
    print np.average(dz)
    dzmin = min(dz)
    
    print("DZ (acc)  : %0.14E"%(w0))
    print("DZ        : %0.14E"%(dzmin))
    
    import matplotlib.pylab as plt
    mdl.plot()
    plt.show()


if __name__ == '__main__':
    #test()
    ex1()
    #ex3()
