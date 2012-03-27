# -*- coding: utf-8 -*-
"""
==========
Kalvomalli
==========

"""

from __future__ import division

import teefem
import teefem.elements
from teefem.boundary_conditions import PressureBoundaryCondition
from common import Model
from teefem import log, cache

#from scipy.sparse.linalg import spsolve
from numpy import array, zeros, matrix
import numpy as np
#import matplotlib.tri as tri
np.set_printoptions(precision = 2, edgeitems=20)

###############################################################################

class MembraneCharacteristic(object):
    ''' Kalvoelementin "ominaissuureita"
    
        Tx,Ty,Txy = esikiristys (prestress)
        funktioina, esim. Tx = lambda k,e: 10
        
        m = massa pinta-alayksikköä kohden
        
        >>> T = lambda k,e: 10
        >>> m1 = lambda k,e: 1
        >>> mcha = MembraneCharestiristic(Tx = T, Ty = T, Txy = T, m = m1)
    
    '''
    def __init__(self,**kwds):
        # Nämä paremmin
        self.Tx = kwds.get('Tx', lambda k,e: 0)
        self.Ty = kwds.get('Ty', lambda k,e: 0)
        self.Txy = kwds.get('Txy', lambda k,e: 0)
        self.m = kwds.get('m', lambda k,e: 0)

memchar = MembraneCharacteristic

###############################################################################

class MEMTR3(teefem.elements.Element2D):
    
    ''' Kolmion muotoinen kalvoelementti '''
    
    def __init__(self, *args, **kwds):
        super(MEMTR3, self).__init__(*args, **kwds)
        
        self.boundary_conditions = set()
        
        self.has_stiffness = True
        self.has_mass = True
        
        self.dimension = 3*1
        self.degrees_of_freedom = ('DZ',)
        self.geom = kwds.get('geom')
#        self.ipoints = [(1.0/3.0, 1.0/3.0)]
#        self.iweights = [(1.0/2.0)]
        self.ipoints = [(0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]
        self.iweights = [(1.0/6.0), (1.0/6.0), (1.0/6.0)]
        
        # Interpoloidaan tuntematonta kenttää samoilla C0-jatkuvilla muoto-
        # funktioilla kuin elementin geometriaa, tässä tapauksessa siis TRIA3
        # muofofunktioilla jotka tulevat suoraan geometrialuokasta
        self.dNdk = self.geom.dNdk
        self.dNde = self.geom.dNde
        self.dNxy = lambda *ke: self.geom.invJ(*ke) * matrix([self.dNdk(*ke),self.dNde(*ke)])
        self.Nxy = self.geom.N
        
    def assign_boundary_condition(self, bc):
        ''' Kalvon painekuorma '''
        if bc.__class__ is PressureBoundaryCondition:
            self.pressure = bc.pressure

    def assign_char(self, char):
        self.char = char
#        if char.__class__ is MembraneCharacteristic:
#            self.Tx = char.Tx
#            self.Ty = char.Ty
#            self.Txy = char.Txy

    @cache
    def Bg(self,*ke):
        ''' Kinemaattinen matriisi '''
        dNxy = self.dNxy(*ke)
        dNdx = dNxy[0]
        dNdy = dNxy[1]
        Bg = matrix(zeros((2,self.dimension)))
        Bg[0,0::1] = dNdx
        Bg[1,0::1] = dNdy
        return Bg

    @cache
    def T(self, *ke):
        ''' "Kalvovoimamatriisi" '''
        Tx = self.char.Tx(*ke)
        Ty = self.char.Ty(*ke)
        Txy = self.char.Txy(*ke)
        return matrix([[Tx,Txy],[Txy,Ty]])

    @property
    def stiffness_matrix(self):
        ''' Jäykkyysmatriisi '''
        K = matrix(zeros((self.dimension,self.dimension)))
        for (W,ke) in zip(self.iweights, self.ipoints):
            detJ = self.detJ(*ke)
            Bg = self.Bg(*ke)
            T = self.T(*ke)
            K += W*Bg.T*T*Bg*detJ
        return K

    @property
    def mass_matrix(self):
        ''' Massamatriisi '''
        M = matrix(zeros((self.dimension,self.dimension)))
        for (W,ke) in zip(self.iweights, self.ipoints):
            detJ = self.detJ(*ke)
            N = matrix(self.geom.N(*ke))
            m = self.char.m(*ke)
            res = W*N.T*m*N*detJ
#            print("N: {0}   m: {1}    detJ: {2}".format(N,m,detJ))
 #           print res
            M += res
        #print M*1000
        return M

    @property
    def force_vector(self):
        ''' Palauttaa kuormitusvektorin '''
        b = array([W*self.pressure(*ke)*self.Nxy(*ke)*self.detJ(*ke) for (W,ke) in zip(self.iweights, self.ipoints)])
        F = zeros(self.dimension)
        F[0::1] = np.sum(b,axis=0)
        return matrix(F).T
        
#    def update(self, U):
#        ''' '''
#        u = zeros(self.dimension)
#        for i in xrange(len(self.nodes)):
#            node = self.nodes[i]
#            node.update_field(
#                field = 'DEPL', 
#                params = self.degrees_of_freedom, 
#                values = (U[node.gdof[0]],))
#            for j in xrange(3):
#                u[i*3+j] = U[node.gdof[j]]
#        self.u = matrix(u).T


###############################################################################

class Membrane(Model):
    
    ''' Kalvomalli. '''
    
    def __init__(self, *args, **kwds):
        super(Membrane, self).__init__(*args, **kwds)
        self.mapping = {
            'Tria3': MEMTR3,
        }
        self.nodedim = 1
        self.nodedofs = ('dz',)
        self.nodeloads = ('fz')
        self.init()

    def plot(self, **kwds):
        ''' Plottaus. '''
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
        
membrane = Membrane

###############################################################################
###############################################################################
###############################################################################

def test1():
    
    ''' Kalvo-ongelma testi, yksikkökalvo jota venytetään toisen akselin suuntaisesti '''

    q = 100
#    w = lambda r: -a**4*q*(1 - r**2/a**2)**2/(64*D) ???

    mesh = teefem.mesh.unitsquare()
    
    mdl = membrane(mesh = mesh)
    
    # Painekuorma
    teefem.assign_bc(
        elements = mdl.elset['OM1'],
        bc = teefem.pressure_bc(pressure = lambda k,e: -q),
        )
    
    # Vasen ja oikea reuna kiinni
    nset1 = set().union(mdl.nset['GA1'], mdl.nset['GA3'])
#    nset1 = set().union(mdl.nset['GA1'], mdl.nset['GA2'], mdl.nset['GA3'], mdl.nset['GA4'])
    teefem.assign_bc(
        nodes = nset1, 
        bc = teefem.dirichlet_bc(encastre = True)
        )

    # Esijännitys
    car = memchar(Tx = lambda k,e: 100)
#    car = memchar(Tx = lambda k,e: 100, Ty = lambda k,e: 100)
    teefem.assign_char(elements = mdl.elements, char = car)

    mdl.static_solve(export_matrices = False)
    
    dz = [node.fields['DEPL']['DZ'] for node in mdl.nodes]
    dymin = min(dz)
    
#    print("DZ (acc)  : %0.14E"%(-w(0)))
    print("DZ        : %0.14E"%(dymin))
#    print("DZ (CA)   : %0.14E"%(-1.74644966293971E-01))
    
    import matplotlib.pylab as plt
    mdl.plot()
    plt.show()


def test2():
    
    ''' Kalvo-ongelma testi, yksikkökalvon ominaismuodot ja taajuudet '''

    # Lähtötiedot
    m = 1
    T = 10
    a = 1
    b = 1

    mesh = teefem.mesh.unitsquare()
    
    mdl = membrane(mesh = mesh)
    
    # Kaikki reunat kiinni
    nset1 = set().union(mdl.nset['GA1'], mdl.nset['GA2'], mdl.nset['GA3'], mdl.nset['GA4'])
    teefem.assign_bc(
        nodes = nset1, 
        bc = teefem.dirichlet_bc(encastre = True)
        )

    # Esijännitys
    car = memchar(Tx = lambda k,e: T, m = lambda k,e: m)
    teefem.assign_char(elements = mdl.elements, char = car)

    mdl.modal_solve(export_matrices = True)

    # Tarkka ratkaisu
    w = lambda i,j: np.sqrt(T/m)*np.sqrt((i*np.pi/a)**2 + (j*np.pi/b)**2)
    print("acc:")
    X = range(1,4)
    for i in X:
        for j in X:
            print("i: {0}  j: {1}  w: {2}".format(i,j,w(i,j)))
#    
#    dz = [node.fields['DEPL']['DZ'] for node in mdl.nodes]
#    dymin = min(dz)
#    
##    print("DZ (acc)  : %0.14E"%(-w(0)))
#    print("DZ        : %0.14E"%(dymin))
##    print("DZ (CA)   : %0.14E"%(-1.74644966293971E-01))
#    
#    import matplotlib.pylab as plt
#    mdl.plot()
#    plt.show()

if __name__ == '__main__':
    #test1()
    test2()
    #ex2()
    #rotsy()
