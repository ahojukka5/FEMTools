# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 02:04:59 2012

@author: Jukka

Mindlin elementti. Tessler ja Hughesin versio.

"""

from __future__ import division

import teefem
from numpy import array, zeros, matrix
import numpy as np
np.set_printoptions(precision = 2, edgeitems=20, threshold=10000, linewidth=150)

import mindlin

class MINTR3TH(mindlin.MINTR3): # Käytetään normaalia Mindlin laattaa pohjana
    
    ''' Mindlin TRIA3 Tessler & Hughes '''
    
    def __init__(self, *args, **kwds):
        super(MINTR3TH, self).__init__(*args, **kwds)
        
        # 'Lisämuotofunktiot'
        self.dLdk = lambda k,e: array([-0.5*e, 0.5*e, 0.0])
        self.dLde = lambda k,e: array([-1.0*e - 0.5*k + 0.5, 0.5*k, 1.0*e - 0.5])
        self.dMdk = lambda k,e: array([-0.5*e - 1.0*k + 0.5, 1.0*k - 0.5, 0.5*e])
        self.dMde = lambda k,e: array([-0.5*k, 0.0, 0.5*k])

    def Bs(self,*ke):
        ''' Kinemaattinen matriisi Bs '''
        dNxy = self.dNxy(*ke)
        dNdx = dNxy[0]
        dNdy = dNxy[1]

        # Lisämuotofunktiot L
        dLdk = self.dLdk(*ke)
        dLde = self.dLde(*ke)
        dLxy = self.invJ(*ke) * matrix([dLdk,dLde])
        dLdx = dLxy[0]
        dLdy = dLxy[1]

        # Lisämuotofunktiot M
        dMdk = self.dMdk(*ke)
        dMde = self.dMde(*ke)
        dMxy = self.invJ(*ke) * matrix([dMdk,dMde])
        dMdx = dMxy[0]
        dMdy = dMxy[1]
        
        Nxy = self.Nxy(*ke)
        Bs = matrix(zeros((2,self.dimension)))
        Bs[0,0::3] = dNdx
        Bs[1,0::3] = dNdy
        Bs[0,2::3] = Nxy
        Bs[1,1::3] = Nxy
        Bs[0,1::3] += dLdx
        Bs[1,1::3] += dLdy
        Bs[0,2::3] += dMdx
        Bs[1,2::3] += dMdy
        return Bs


class MINTH(mindlin.MindlinModel):
    ''' Reissner-Mindlin Tessler & Hughes model '''
    def __init__(self, *args, **kwds):
        super(MINTH, self).__init__(*args, **kwds)
        self.mapping = {
            'Tria3': MINTR3TH,
        }
        self.init()
    

minth = MINTH

###############################################################################
###############################################################################
###############################################################################

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
    mdl = MINTH(mesh = mesh)
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


def genshapef():
    ''' Tessler ja Hughes muotofunktioiden pyörittelyä '''
    import sympy
    k,e = sympy.var('k,e')
    
    x1 = 0
    y1 = 0
    x2 = 1
    y2 = 0
    x3 = 0
    y3 = 1
    
    a1 = x2-x3
    a2 = x3-x1
    a3 = x1-x2
    b1 = y3-y2
    b2 = y1-y3
    b3 = y2-y1
    
    print a1,a2,a3,b1,b2,b3
    
    N1 = 1-k-e
    N2 = k
    N3 = e
    
    NN1 = N1*(2*N1-1)
    NN2 = N2*(2*N2-1)
    NN3 = N3*(2*N3-1)
    NN4 = 4*N1*N2
    NN5 = 4*N2*N3
    NN6 = 4*N3*N1
    
    print NN1.expand()
    print NN2.expand()
    print NN3.expand()
    print NN4.expand()
    print NN5.expand()
    print NN6.expand()
    
    L1 = 1/8*(b3*NN4 - b2*NN6)
    L2 = 1/8*(b1*NN5 - b3*NN4)
    L3 = 1/8*(b2*NN6 - b1*NN5)
    M1 = 1/8*(a2*NN6 - a3*NN4)
    M2 = 1/8*(a3*NN4 - a1*NN5)
    M3 = 1/8*(a1*NN5 - a2*NN6)
    
    dL1dk = L1.diff(k,1)
    dL2dk = L2.diff(k,1)
    dL3dk = L3.diff(k,1)
    dL1de = L1.diff(e,1)
    dL2de = L2.diff(e,1)
    dL3de = L3.diff(e,1)
    
    dM1dk = M1.diff(k,1)
    dM2dk = M2.diff(k,1)
    dM3dk = M3.diff(k,1)
    dM1de = M1.diff(e,1)
    dM2de = M2.diff(e,1)
    dM3de = M3.diff(e,1)
    
    dLdk = sympy.Matrix([dL1dk, dL2dk, dL3dk])
    dLde = sympy.Matrix([dL1de, dL2de, dL3de])
    dMdk = sympy.Matrix([dM1dk, dM2dk, dM3dk])
    dMde = sympy.Matrix([dM1de, dM2de, dM3de])
    
    print("dLdk:\n{0}".format(dLdk))
    print("dLde:\n{0}".format(dLde))
    print("dMdk:\n{0}".format(dMdk))
    print("dMde:\n{0}".format(dMde))

if __name__ == '__main__':
    ex1()
    #ex2()
    #ex3()
