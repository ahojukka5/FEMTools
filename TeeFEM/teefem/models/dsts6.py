# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 02:04:59 2012

@author: Jukka

Mindlin elementti, perusversio.

"""

from __future__ import division

import teefem
import numpy as np
from numpy import array, zeros, matrix
import common

import mindlin

class DSTS6TR3(mindlin.MINTR3): # Lähdetään liikkeelle perus Mindlin elementistä
    
    ''' DSTS6TR3 TRIA3 Reissner-Mindlin finite element '''
    
    def __init__(self, *args, **kwds):
        super(DSTS6TR3, self).__init__(*args, **kwds)
        self.dimension = 3*5
        self.degrees_of_freedom = ('DZ','DRX','DRY','GX','GY')

        # Kolmen pisteen integroimiskaava
#        self.ipoints = [(0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]
#        self.iweights = [(1.0/6.0), (1.0/6.0), (1.0/6.0)]

        # Kolmen pisteen integroimiskaava (Hammer)
#        k = 1.0/6.0
#        self.ipoints = [(k, k), (2*k, k), (k, 2*k)]
#        self.iweights = [(0.5*k), (0.5*k), (0.5*k)]

        self.ipoints = [(1.0/3.0, 1.0/3.0)]
        self.iweights = [(1.0/2.0)]

        # Muotofunktiot
        self.Lxy = lambda k,e: array([1-k-e,k,e])

        dRdk = lambda k,e: array([6.0*e, -3.0*e, -3.0*e])
        dRde = lambda k,e: array([12.0*e + 6.0*k - 6.0, -3.0*k, -12.0*e - 3.0*k + 6.0])
        dRxdk = lambda k,e: array([3.0*e - 1.0, -1.5*e + 1.0, 1.5*e])
        dRxde = lambda k,e: array([6.0*e + 3.0*k - 4.0, -1.5*k, 6.0*e + 1.5*k - 2.0])
        dRydk = lambda k,e: array([0, -1.5*e, -1.5*e])
        dRyde = lambda k,e: array([0, -1.5*k, -1.5*k])
        dRgdk = lambda k,e: array([3.0*e, -1.5*e, 1.5*e])
        dRgde = lambda k,e: array([6.0*e + 3.0*k - 3.0, -1.5*k, 6.0*e + 1.5*k - 3.0])
        dHdk = lambda k,e: array([6.0*e + 12.0*k - 6.0, -9.0*e - 12.0*k + 6.0, 3.0*e])
        dHde = lambda k,e: array([6.0*k, -9.0*k, 3.0*k])
        dHxdk = lambda k,e: array([0, -1.5*e, -1.5*e])
        dHxde = lambda k,e: array([0, -1.5*k, -1.5*k])
        dHydk = lambda k,e: array([3.0*e + 6.0*k - 4.0, 1.5*e + 6.0*k - 2.0, -1.5*e])
        dHyde = lambda k,e: array([3.0*k - 1.0, 1.5*k, -1.5*k + 1.0])
        dHgdk = lambda k,e: array([3.0*e + 6.0*k - 3.0, 1.5*e + 6.0*k - 3.0, -1.5*e])
        dHgde = lambda k,e: array([3.0*k, 1.5*k, -1.5*k])

        invJ = self.geom.invJ
        self.dRxy = lambda *ke: invJ(*ke) * matrix([dRdk(*ke),dRde(*ke)])
        self.dRxxy = lambda *ke: invJ(*ke) * matrix([dRxdk(*ke),dRxde(*ke)])
        self.dRyxy = lambda *ke: invJ(*ke) * matrix([dRydk(*ke),dRyde(*ke)])
        self.dRgxy = lambda *ke: invJ(*ke) * matrix([dRgdk(*ke),dRgde(*ke)])
        self.dHxy = lambda *ke: invJ(*ke) * matrix([dHdk(*ke),dHde(*ke)])
        self.dHxxy = lambda *ke: invJ(*ke) * matrix([dHxdk(*ke),dHxde(*ke)])
        self.dHyxy = lambda *ke: invJ(*ke) * matrix([dHydk(*ke),dHyde(*ke)])
        self.dHgxy = lambda *ke: invJ(*ke) * matrix([dHgdk(*ke),dHgde(*ke)])
        self.detJ = self.geom.detJ


#    def Db(self,*ke):
#        E = self.material.E
#        nu = self.material.nu
#        return E*self.thickness(*ke)**3/(12*(1-nu**2)) * matrix([[1,nu,0],[nu,1,0],[0,0,0.5*(1-nu)]])
#
#    def Ds(self,*ke):
#        G = self.material.E
#        nu = self.material.nu
#        return 5*E*self.thickness(*ke)/(12*(1+nu)) * matrix([[1,0],[0,1]])
#

    def Bb(self,*ke):
        ''' Kinemaattinen matriisi Bb '''
        dR = self.dRxy(*ke)
        dRx = self.dRxxy(*ke)
        dRy = self.dRyxy(*ke)
        dRg = self.dRgxy(*ke)
        dH = self.dHxy(*ke)
        dHx = self.dHxxy(*ke)
        dHy = self.dHyxy(*ke)
        dHg = self.dHgxy(*ke)
        Bb = matrix(zeros((3,self.dimension)))
        x = 0
        y = 1
        Bb[0,0::5] = dR[x]
        Bb[0,1::5] = dRx[x]
        Bb[0,2::5] = dRy[x]
        Bb[0,3::5] = dRg[x]
        Bb[0,4::5] = dRy[x]
        Bb[1,0::5] = dH[y]
        Bb[1,1::5] = dHx[y]
        Bb[1,2::5] = dHy[y]
        Bb[1,3::5] = dHx[y]
        Bb[1,4::5] = dHg[y]
        Bb[2,0::5] = dR[y]+dH[x]
        Bb[2,1::5] = dRx[y]+dHx[x]
        Bb[2,2::5] = dRy[y]+dHy[x]
        Bb[2,3::5] = dRg[y]+dHx[x]
        Bb[2,4::5] = dRy[y]+dHg[x]
        return Bb

    def Bs(self,*ke):
        ''' Kinemaattinen matriisi Bs '''
        Lxy = self.Lxy(*ke)
        Bs = matrix(zeros((2,self.dimension)))
        Bs[0,3::5] = Lxy
        Bs[1,4::5] = Lxy
        return Bs

#    @property
#    def stiffness_matrix(self):
#        ''' Jäykkyysmatriisi '''
#        Kb = matrix(zeros((self.dimension,self.dimension)))
#        Ks = matrix(zeros(Kb.shape))
#        for (W,ke) in zip(self.iweights, self.ipoints):
#            detJ = self.detJ(*ke)
#            Bb = self.Bb(*ke)
#            Db = self.Db(*ke)
#            Kb += W*Bb.T*Db*Bb*detJ
#            Bs = self.Bs(*ke)
#            Ds = self.Ds(*ke)
#            Ks += W*Bs.T*Ds*Bs*detJ
#        return Kb + Ks

    @property
    def force_vector(self):
        ''' Palauttaa kuormitusvektorin '''
        b = array([W*self.pressure(*ke)*self.Lxy(*ke)*self.detJ(*ke) for (W,ke) in zip(self.iweights, self.ipoints)])
        F = zeros(self.dimension)
        F[0::5] = np.sum(b,axis=0)
        return matrix(F).T

#    def update(self, U):
#        u = zeros(self.dimension)
#        for i in xrange(len(self.nodes)):
#            node = self.nodes[i]
#            node.update_field(
#                field = 'DEPL', 
#                params = self.degrees_of_freedom, 
#                values = (U[node.gdof[0]], U[node.gdof[1]], U[node.gdof[2]]))
#            for j in xrange(3):
#                u[i*3+j] = U[node.gdof[j]]
#        self.u = matrix(u).T

class MINDSTS6(common.Model):
    ''' DST-S6 '''

    def __init__(self, *args, **kwds):
        super(MINDSTS6, self).__init__(*args, **kwds)
        self.nodedim = 5
        self.nodedofs = ('dz','drx','dry','gx','gy')
        self.nodeloads = ('fz','mx','my','gx','gy')
        self.mapping = {
            'Tria3': DSTS6TR3,
        }
        self.init()


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

dsts6 = MINDSTS6

###############################################################################
###############################################################################
###############################################################################



def test():
    
    ''' DSTS6TR3-elementin matriisien tarkistus '''
    
    n1 = teefem.geom.Node(x=0, y=0, z=0.0)
    n2 = teefem.geom.Node(x=1, y=0, z=0.0)
    n3 = teefem.geom.Node(x=0, y=1, z=0.0)

    g1 = teefem.geom.Tria3(nodes = (n1,n2,n3))
    g1.material = teefem.materials.Elastic(E = 100.0e9, nu = 0.3)

    e1 = DSTS6TR3(geom = g1)
    e1.pressure = lambda *ke: 100e3
    e1.thickness = lambda *ke: 10e-3

    print e1.status

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
    mdl = MINDSTS6(mesh = mesh)
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


def shapef():
    ''' Pyörittelee DST-6 funktiot '''
    import sympy
    k,e = sympy.var('k,e')
    x1 = 0
    y1 = 0
    x2 = 1
    y2 = 0
    x3 = 0
    y3 = 1
    c1 = x2*y3 - x3*y2
    c2 = x3*y1 - x1*y3
    c3 = x1*y2 - x2*y1
    b1 = y2-y3
    b2 = y3-y1
    b3 = y1-y2
    a1 = x3-x2
    a2 = x1-x3
    a3 = x2-x1
    x21 = x2-x1
    y31 = y3-y1
    y21 = y2-y1
    x31 = x3-x1
    A2 = x21*y31 - y21*x31
    L1 = (c1+b1*k+a1*e)/A2
    L2 = (c2+b2*k+a2*e)/A2
    L3 = (c3+b3*k+a3*e)/A2
    print a1,a2,a3
    print b1,b2,b3
    print c1,c2,c3
    print("L1: {0}".format(L1))
    print("L2: {0}".format(L2))
    print("L3: {0}".format(L3))
    N1 = (2*L1-1)*L1
    N2 = (2*L2-1)*L2
    N3 = (2*L3-1)*L3
    N4 = 4*L1*L2
    N5 = 4*L2*L3
    N6 = 4*L3*L1
    print("N1: {0}".format(N1))
    print("N2: {0}".format(N2))
    print("N3: {0}".format(N3))
    print("N4: {0}".format(N4))
    print("N5: {0}".format(N5))
    print("N6: {0}".format(N6))
    x21 = x2-x1
    y21 = y2-y1
    x32 = x3-x2
    y32 = y3-y2
    x13 = x1-x3
    y13 = y1-y3
    l4 = sympy.sqrt(x21**2 + y21**2)
    l5 = sympy.sqrt(x32**2 + y32**2)
    l6 = sympy.sqrt(x13**2 + y13**2)
    print("l4: {0}".format(l4))
    print("l5: {0}".format(l5))
    print("l6: {0}".format(l6))
    # ???
#    m4 = x21/l4
#    n4 = y21/l4
#    m5 = x32/l5
#    n5 = y32/l5
#    m6 = x13/l6
#    n6 = y13/l6
    m4 = 0
    n4 = 1
    m5 = sympy.sqrt(2)/2
    n5 = sympy.sqrt(2)/2
    m6 = -1
    n6 = 0
    print("4: {0}".format((m4,n4)))
    print("5: {0}".format((m5,n5)))
    print("6: {0}".format((m6,n6)))
    R1 = 1.5*(m6*N6/l6 - m4*N4/l4)
    R2 = 1.5*(m4*N4/l4 - m5*N5/l5)
    R3 = 1.5*(m5*N5/l5 - m6*N6*l6)
    print("R1: {0}".format(R1))
    print("R2: {0}".format(R2))
    print("R3: {0}".format(R3))
    Rx1 = N1 + N4*(0.5*n4**2 - 0.25*m4**2) + N6*(0.5*n6**2 - 0.25*m6**2)
    Rx2 = N2 + N4*(0.5*n4**2 - 0.25*m4**2) + N5*(0.5*n5**2 - 0.25*m5**2)
    Rx3 = N3 + N5*(0.5*n5**2 - 0.25*m5**2) + N6*(0.5*n6**2 - 0.25*m6**2)
    print("Rx1: {0}".format(Rx1))
    print("Rx2: {0}".format(Rx2))
    print("Rx3: {0}".format(Rx3))
    Ry1 = -0.75*(m4*n4*N4 + m6*n6*N6)
    Ry2 = -0.75*(m4*n4*N4 + m5*n5*N5)
    Ry3 = -0.75*(m5*n5*N5 + m6*n6*N6)
    print("Ry1: {0}".format(Ry1))
    print("Ry1: {0}".format(Ry2))
    print("Ry1: {0}".format(Ry3))
    H1 = 1.5*(n6*N6/l6 - n4*N4/l4)
    H2 = 1.5*(n4*N4/l4 - n5*N5/l5)
    H3 = 1.5*(n5*N5/l5 - n6*N6/l6)
    print("H1: {0}".format(H1))
    print("H2: {0}".format(H2))
    print("H3: {0}".format(H3))
    Hx1 = Ry1
    Hx2 = Ry2
    Hx3 = Ry3
    print("Hx1: {0}".format(Hx1))
    print("Hx2: {0}".format(Hx2))
    print("Hx3: {0}".format(Hx3))
    Hy1 = N1 + N4*(0.5*m4**2 - 0.25*n4**2) + N6*(0.5*m6**2 - 0.25*n6**2)
    Hy2 = N2 + N4*(0.5*m4**2 - 0.25*n4**2) + N5*(0.5*m5**2 - 0.25*n5**2)
    Hy3 = N3 + N5*(0.5*m5**2 - 0.25*n5**2) + N6*(0.5*m6**2 - 0.25*n6**2)
    print("Hy1: {0}".format(Hy1))
    print("Hy2: {0}".format(Hy2))
    print("Hy3: {0}".format(Hy3))
    Rg1 = Rx1 - L1
    Hg1 = Hy1 - L1
    print("Rg1: {0}".format(Rg1))
    print("Hg1: {0}".format(Hg1))
    Rg2 = Rx2 - L2
    Hg2 = Hy2 - L2
    print("Rg2: {0}".format(Rg2))
    print("Hg2: {0}".format(Hg2))
    Rg3 = Rx3 - L3
    Hg3 = Hy3 - L3
    print("Rg3: {0}".format(Rg3))
    print("Hg3: {0}".format(Hg3))
    R = [R1,R2,R3]
    Rx = [Rx1,Rx2,Rx3]
    Ry = [Ry1,Ry2,Ry3]
    Rg = [Rg1,Rg2,Rg3]
    H = [H1,H2,H3]
    Hx = [Hx1,Hx2,Hx3]
    Hy = [Hy1,Hy2,Hy3]
    Hg = [Hg1,Hg2,Hg3]
    print("dRdk = lambda k,e: array({0})".format([Ri.diff(k) for Ri in R]))
    print("dRde = lambda k,e: array({0})".format([Ri.diff(e) for Ri in R]))
    print("dRxdk = lambda k,e: array({0})".format([Ri.diff(k) for Ri in Rx]))
    print("dRxde = lambda k,e: array({0})".format([Ri.diff(e) for Ri in Rx]))
    print("dRydk = lambda k,e: array({0})".format([Ri.diff(k) for Ri in Ry]))
    print("dRyde = lambda k,e: array({0})".format([Ri.diff(e) for Ri in Ry]))
    print("dRgdk = lambda k,e: array({0})".format([Ri.diff(k) for Ri in Rg]))
    print("dRgde = lambda k,e: array({0})".format([Ri.diff(e) for Ri in Rg]))
    print("dHdk = lambda k,e: array({0})".format([Hi.diff(k) for Hi in H]))
    print("dHde = lambda k,e: array({0})".format([Hi.diff(e) for Hi in H]))
    print("dHxdk = lambda k,e: array({0})".format([Hi.diff(k) for Hi in Hx]))
    print("dHxde = lambda k,e: array({0})".format([Hi.diff(e) for Hi in Hx]))
    print("dHydk = lambda k,e: array({0})".format([Hi.diff(k) for Hi in Hy]))
    print("dHyde = lambda k,e: array({0})".format([Hi.diff(e) for Hi in Hy]))
    print("dHgdk = lambda k,e: array({0})".format([Hi.diff(k) for Hi in Hg]))
    print("dHgde = lambda k,e: array({0})".format([Hi.diff(e) for Hi in Hg]))

def shapeftest():
    def tsk(a): # Tangenttivektorin i,j yksikkövektorit:
        i = np.cos(a)
        j = np.sin(a)
        return "theta sk kulma {0} = {1}i + {2}j".format(a/np.pi*180,i,j)
    def tnk(a): # Normaalivektorin i,j yksikkövektorit:
        i = -np.sin(a)
        j = np.cos(a)
        return "theta nk kulma {0} = {1}i + {2}j".format(a/np.pi*180,i,j)
    # Tangentit
    print tsk(0) # Pitäisi olla 1i + 0j
    print tsk(3/4*np.pi) # Pitäisi olla -sqrt(2)/2 i + sqrt(2)/2 j
    print tsk(-np.pi/2) # Pitäisi olla 0i - 1j
    # Normaalit
    print tnk(0) # Pitäisi olla 0i -1j
    print tnk(3/4*np.pi) # Pitäisi olla sqrt(2)/2 i + sqrt(2)/2 j
    print tnk(-np.pi/2) # Pitäisi olla -1i + 0j

if __name__ == '__main__':
    #shapef()
    #shapeftest()
    #test()
    ex1()
    #ex3()
