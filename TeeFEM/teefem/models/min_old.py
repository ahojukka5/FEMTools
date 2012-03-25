# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 02:04:59 2012

@author: Jukka

Mindlin elementti, perusversio.

Vanha versio.

"""

from __future__ import division
from scipy.sparse.linalg import spsolve
import teefem
import logging
from numpy import array, zeros, concatenate, linspace, matrix
from scipy.sparse.linalg import spsolve
from numpy import sqrt
from common import Model
import numpy as np
import matplotlib.tri as tri
np.set_printoptions(precision = 2, edgeitems=20, threshold=10000, linewidth=150)

def checkmatrix(A):
    teefem.log.info('Determinant: {0}'.format(np.linalg.det(A)))
    s = A.shape
    teefem.log.info('Shape: {0}'.format(s))
    for i in range(s[0]):
        print("Sum col {0}: {1}".format(i+1, np.sum(A[i,:])))
        print("Sum row {0}: {1}".format(i+1, np.sum(A[:,i])))

def readmailfile(fh):
    ''' Reads mesh (stream) in format Aster mail '''

def encastre(*args):
    ''' Encastre boundary condition. '''
    return 0

class Force_coque(object):
    ''' FORCE_COQUE boundary condition. Uniform pressure to plate. '''
    def __init__(self, *args, **kwds):
        self.elements = kwds.get('group_ma')
        self.pressure = kwds.get('pressure')

class DDL_impo(object):
    ''' DDL_IMPO boundary condition.'''
    def __init__(self, *args, **kwds):
        self.nodes = kwds.get('group_no')
        self.pressure = kwds.get('pressure')
        self.liaison = kwds.get('liaison')

class MINTR3(teefem.elements.Element2D):
    ''' MINTR3 TRIA3 Reissner-Mindlin finite element '''
    def __init__(self, *args, **kwds):
        super(MINTR3, self).__init__(*args, **kwds)
        self.dimension = 3*3
        self.has_stiffness = True

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

        # Nämä paremmin
        self.thickness = lambda *ke: 0.500
        self.kappa = 5/6
        self.pressure = lambda *ke: 100e3

        # Materiaaliominaisuudet. Paksuus on funktio.
        self.material = self.geom.material
        E = self.material.E
        G = self.material.G
        nu = self.material.nu
        self.Db = lambda *ke: E*self.thickness(*ke)**3/(12*(1-nu**2)) * matrix([[1,nu,0],[nu,1,0],[0,0,0.5*(1-nu)]])
        self.Ds = lambda *ke: self.kappa * G * self.thickness(*ke) * matrix([[1,0],[0,1]])

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

class MINTR3SI(MINTR3):
    ''' MINTR3 TRIA3 Reissner-Mindlin finite element, reduced integration '''
    def __init__(self, *args, **kwds):
        super(MINTR3SI, self).__init__(*args, **kwds)
        self.ipoints = [(1.0/3.0, 1.0/3.0)]
        self.iweights = [(1.0/2.0)]

class MIN3(Model):
    
    ''' Reissner-Mindlin model, reduced integration '''
    
    def __init__(self, *args, **kwds):
        super(MIN3, self).__init__(*args, **kwds)
        mapping = {
#            'Seg2': MEBODST,
            'Tria3': MINTR3,
        }
        # Create elements according element mapping
        for (name, net) in self.mesh.nets.iteritems():
            try:
                shape = net.__class__.__name__
                self.elements[name] = mapping[shape](geom = net)
            except KeyError:
                logging.debug("No corresponding finite element to net {0} in DKT model".format(shape))
        self.group_ma = {}
        for (k,v) in self.mesh.group_ma.iteritems():
            try: 
                self.group_ma[k] = set([self.elements[net.name] for net in v])
                teefem.log.info("New finite element group: {0}".format(k))
            except KeyError:
                teefem.log.warning("Unable to create finite element group for group_ma {0}".format(k))
        teefem.log.info("Model initialized. Number of elements: {0}".format(len(self.elements)))

    def assembly_stiffness_matrix(self):
        ''' Return DKT global stiffness matrix '''
        
        # 0) Size of stiffness matrix
        ijvsize = 0
        kdim = 0
        nodedim = 3
        teefem.log.debug("Number of elements: {0}".format(len(self.elements)))
        for element in self.elements.itervalues():
            if not element.has_stiffness:
                continue
            ijvsize += element.dimension**2
            kdim += element.dimension
            
        # 1) NUME_DDL (assign position to global stiffness matrix)
        idx = 0
        for node in self.nodes.itervalues():
            node.gdof = xrange(idx*nodedim,(idx+1)*nodedim)
            if hasattr(node, 'bc'):
                if node.bc.has_key('dx'):
                    ijvsize += 2
                if node.bc.has_key('dy'):
                    ijvsize += 2
                if node.bc.has_key('liaison'):
                    ijvsize += 3*2
            idx += 1
        ijv = teefem.common.IJV(ijvsize)
        
#        for node in self.nodes.itervalues():
#            print node.gdof
        
        # 2) Assembly element local stiffness matrix to global stiffness matrix
        for element in self.elements.itervalues():
            if not element.has_stiffness:
                continue
            k_loc = element.stiffness_matrix
#            print k_loc
            gdof = array([node.gdof for node in element.nodes]).flatten()
#            print gdof
            for I in xrange(element.dimension):
                for J in xrange(element.dimension):
                    ijv.add(gdof[I], gdof[J], k_loc[I,J])
                    
        # 3) Assign boundary conditions to stiffness matrix with Lagrange multipliers
        alp = 1e6
        lidx = len(self.nodes)*3
        print("lidx: {0}".format(lidx))
#        for node in self.nodes.itervalues():
#            if hasattr(node, 'bc'):
#                if node.bc.has_key('dx'):
#                    ijv.add(node.gdof[0],lidx,alp)
#                    ijv.add(lidx,node.gdof[0],alp)
#                    lidx += 1
#                if node.bc.has_key('dy'):
#                    ijv.add(node.gdof[1],lidx,alp)
#                    ijv.add(lidx,node.gdof[1],alp)
#                    lidx += 1

        # Dirty hack
        for node in self.nodes.itervalues():
            if hasattr(node, 'bc'):
                if node.bc.has_key('liaison'):
                    f = np.sqrt(node.x**2 + node.y**2 + node.z**2)
                    teefem.log.info('Locking node ({0},{1},{2}, r={3})'.format(node.x, node.y, node.z, f))
                    ijv.add(node.gdof[0],lidx,alp)
                    ijv.add(lidx,node.gdof[0],alp)
                    lidx += 1
                    ijv.add(node.gdof[1],lidx,alp)
                    ijv.add(lidx,node.gdof[1],alp)
                    lidx += 1
                    ijv.add(node.gdof[2],lidx,alp)
                    ijv.add(lidx,node.gdof[2],alp)
                    lidx += 1
                    
        self.total_dimension = lidx # Total dimension of stiffness matrix
        return ijv.tocoo()
    
    def assembly_force_vector(self):
        ''' Return DKT force vector '''
        R = zeros(self.total_dimension)

        for element in self.elements.itervalues():
            if not hasattr(element, 'bc'):
                continue
            bc = element.bc
            element.qn = lambda k: bc['pressure']
            element.qt = lambda k: bc['shear']
            fq = element.force_vector
            gdof = array([node.gdof for node in element.nodes]).flatten()
            for I in xrange(element.dimension):
                R[gdof[I]] += fq[I]
        lidx = len(self.nodes)*3

        for node in self.nodes.itervalues():
            if hasattr(node, 'bc'):
                if node.bc.has_key('fx'):
                    R[node.gdof[0]] += node.bc['fx']
                if node.bc.has_key('fy'):
                    R[node.gdof[1]] += node.bc['fy']
                if node.bc.has_key('dx'):
                    R[lidx] = node.bc['dx']*1e12
                    lidx += 1
                if node.bc.has_key('dy'):
                    R[lidx] = node.bc['dy']*1e12
                    lidx += 1
        return R
        
    def static_solve(self, **kwds):
        debug = kwds.get('print_matrices', False)
        teefem.log.info("Creating stiffness matrix")
        K = self.assembly_stiffness_matrix().tocsr()
#        if debug: print K.todense()
        teefem.log.info("Creating force vector")
        R = self.assembly_force_vector()
#        if debug: print R
        teefem.log.info("Solving")
#        checkmatrix(K.todense())
        U = spsolve(K,R)
#        ndim = len(self.nodes)*3
#        print("Max DZ: {0}".format(max(U[0:ndim:3])))
#        if debug: 
#            print U
        if kwds.get('print_matrices', False): 
            np.savetxt('K.txt', K.todense(), fmt='%+05.2E')
            np.savetxt('R.txt', R, fmt='%+05.2E')
            np.savetxt('U.txt', R, fmt='%+05.2E')
        teefem.log.info("Updating fields")
        for node in self.nodes.itervalues():
            node.update_field(
                field = 'DEPL', 
                params = ('DZ','DRX','DRY'), 
                values = (U[node.gdof[0]], U[node.gdof[1]], U[node.gdof[2]]),
                )

###############################################################################
###############################################################################
###############################################################################



def test():
    
    ''' MINTR3SI-elementin matriisien tarkistus '''
    
    n1 = teefem.geom.Node(x=0, y=0, z=0.0)
    n2 = teefem.geom.Node(x=1, y=0, z=0.0)
    n3 = teefem.geom.Node(x=0, y=1, z=0.0)

    g1 = teefem.geom.Tria3(nodes = (n1,n2,n3))
    g1.material = teefem.materials.Elastic(E = 100.0e9, nu = 0.3)

    e1 = MINTR3SI(geom = g1)
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

def ex3():
    
    ''' DKT esimerkki 3. Ratkaistaan suorakaidelaattaa, neljä DKT-elementtiä
    
        4_____3
        |\3 /|
        | \/ |
        |4/\2| (solmupiste 5 keskellä)
        |/1_\|
        1    2
    '''

    L = 1
    H = 1    
    t = 10e-3
    q = 100e3

    n1 = teefem.geom.Node(x=0, y=0, z=0.0)
    n2 = teefem.geom.Node(x=L, y=0, z=0.0)
    n3 = teefem.geom.Node(x=L, y=H, z=0.0)
    n4 = teefem.geom.Node(x=0, y=H, z=0.0)
    n5 = teefem.geom.Node(x=L/2, y=H/2, z=0.0)
    g1 = teefem.geom.Tria3(name = 'M1', nodes = (n1,n2,n5))
    g2 = teefem.geom.Tria3(name = 'M2', nodes = (n2,n3,n5))
    g3 = teefem.geom.Tria3(name = 'M3', nodes = (n3,n4,n5))
    g4 = teefem.geom.Tria3(name = 'M4', nodes = (n4,n1,n5))
    mesh = teefem.geom.Mesh()
    mesh.group_ma['OM1'] = set([g1,g2,g3,g4])
    mesh.group_no['GA1'] = set([n1,n2,n3,n4])
    print mesh.group_ma
    print mesh.group_no
    mesh.nodes['N1'] = n1
    mesh.nodes['N2'] = n2
    mesh.nodes['N3'] = n3
    mesh.nodes['N4'] = n4
    mesh.nodes['N5'] = n5
    mesh.nets['M1'] = g1
    mesh.nets['M2'] = g2
    mesh.nets['M3'] = g3
    mesh.nets['M4'] = g4

    mdl = MEDSTR3(mesh = mesh)
    print mdl.group_ma
    print mdl.group_no
    mat = teefem.materials.Elastic(E = 210.0e9, nu = 0.3)
    mesh.assign_material(mat)

    mdl.assign_bc(elem_group = 'OM1', pressure = q)
    mdl.assign_bc(node_group = 'GA1', liaison='encastre')

    mdl.static_solve(print_matrices = True)
    
    dymin = 0
    
    for node in mdl.nodes.itervalues():
        dz = node.fields['DEPL']['DZ']
        print dz
        if dz > dymin: dymin = dz
    
    print("DZ     : %0.14E"%(dymin))

def ex1():
    
    ''' Mindlin esimerkki 1, pyörähdyssymmetrinen laatta jakautuneella kuormalla '''
    
    meshfile = 'tria3.mail' # Linear triangle mesh file
    
    t = 500e-3
    q = 100e3
    nu = 0.3
    E = 100.0e9
    a = 1
    D = E*t**3/(12*(1-nu**2))
    phi = 16/5*(t/a)**2/(1-nu)
    print phi
    w0 = q*a**4/(64*D)*(1+phi)

    mesh = None
    
    with open(meshfile,'r') as fh:
        mesh = teefem.geom.parsemail(fh)

    mat = teefem.materials.Elastic(E = E, nu = nu)
    mesh.assign_material(mat)
    
    mdl = MIN3(mesh = mesh)
    
#    om1 = mesh.find_group_ma('OM1')    
#    bc1 = Force_coque(group_ma = om1, pressure = q)
#    mdl.assign_bc(bc1)
#
#    ga1 = mesh.find_group_no('GA1')
#    bc2 = DDL_impo(group_no = ga1, liaison = encastre)
#    mdl.assign_bc(bc2)

    mdl.assign_bc(elem_group = 'OM1', pressure = -q)
    mdl.assign_bc(node_group = 'GA1', liaison='encastre')

    mdl.static_solve(print_matrices = False)
   
    dymin = 0
    
    for node in mdl.nodes.itervalues():
        dz = node.fields['DEPL']['DZ']
#        print dz
        if dz > dymin: dymin = dz
    
    print("DZ (TF)  : %0.14E"%(dymin))
#    print("DZ (CA)  : %0.14E"%(1.74723941631249E-01)) # 10 mm
    print("DZ (CA)  : %0.14E"%(2.99736969818759E-06)) # 500 mm
    print("DZ (acc) : %0.14E"%(w0))

    import matplotlib.pylab as plt
    mdl.plot()
    plt.show()

if __name__ == '__main__':
    #test()
    ex1()
    #ex3()
