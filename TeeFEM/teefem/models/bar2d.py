# -*- coding: utf-8 -*-
"""
    2d sauvamalli
    
    kesken.
    
"""
from __future__ import division
import teefem
import logging
import numpy as np
from common import Model
from numpy import array, zeros, matrix
from numpy import sqrt

###############################################################################

#class Bar(teefem.elements.Element12):
#    '''
#    Plane stress 1d elements, general class.
#    '''
#    
#    def __init__(self, *args, **kwds):
#        super(Bar, self).__init__(*args, **kwds)
#
#    @property
#    def force_vector(self):
#        Jk = self.geom.Jk
#        Je = self.geom.Je
#        N = self.geom.N
#        Rq = zeros(self.dimension)
#        qx = array([W * self.qn(k) * N(k) * Jk(k) for (W,k) in zip(self.iweights,self.ipoints)])
#        qy = array([W * self.qn(k) * N(k) * Je(k) for (W,k) in zip(self.iweights,self.ipoints)])
#        Rq[1::2] += np.sum(qy,axis=0)
#        Rq[0::2] += np.sum(qx,axis=0)
#        return Rq


class Bar2D(object):
    ''' Sauvaelementti '''
    def __init__(self, **kwds):
        self.nodes = kwds['nodes']
        self.ipoints = [(0.0)]
        self.iweights = [(2.0)]
        self.dimension = 4
        self.has_stiffness = True
        self.boundary_conditions = set()
        self.degrees_of_freedom = ('DX','DY')
        self.geom = kwds.get('geom')
        self.dNdk = self.geom.dNdk
        self.dNde = self.geom.dNde

    def assign_material(self, mat):
        self.material = mat
#        log.debug("Material assigned. E: {0}".format(self.material.E))

#    def assign_boundary_condition(self, bc):
#        if bc.__class__ is PressureBoundaryCondition:
#            self.pressure = bc.pressure

    def assign_char(self, char):
        if char.__class__ is PlateCharacteristic:
            self.thickness = char.thickness
        
    @property
    def stiffness_matrix(self):
        alp = teefem.anglxy(self.nodes[0], self.nodes[1])
        he = teefem.dist(self.nodes[0], self.nodes[1])
        c = np.cos(alp)
        s = np.sin(alp)
        cc = c*c
        cs = c*s
        ss = s*s
        return self.E*self.A/he*np.matrix([[cc,cs,-cc,-cs],[cs,ss,-cs,-ss],[-cc,-cs,cc,cs],[-cs,-ss,cs,ss]])

    @property
    def detJ(self):
        return self.geom.detJ

    def B(self,*ke):
        """
        Returns plane kinematic matrix B 
        """
        pass

#    def update(self, U):
#        ''' '''
#        pass
#        u = zeros(self.dimension)
#        for i in xrange(len(self.nodes)):
#            node = self.nodes[i]
#            node.update_field(
#                field = 'DEPL', 
#                params = ('DX','DY'), 
#                values = (U[node.gdof[0]], U[node.gdof[1]]))
#            for j in xrange(2):
#                u[i*2+j] = U[node.gdof[j]]
#        self.u = matrix(u).T

###############################################################################

class Bar2DModel(Model):
    
    ''' Bar 2d model '''
    
    def __init__(self, *args, **kwds):
        super(Bar2DModel, self).__init__(*args, **kwds)
        self.mapping = {
            'Seg2': Bar2D,
        }
        self.nodedim = 2
        self.nodedofs = ('dx','dy')
        self.nodeloads = ('fz','fy')
        self.init()
    


logging.debug("Module {0} loaded.".format(__file__))

###############################################################################

def ex1():
    ''' Ristikkorakenteen ratkaisu  '''
    
    n1 = teefem.geom.Node(x=0, y=0, z=0)
    n2 = teefem.geom.Node(x=1, y=0, z=0)
    n3 = teefem.geom.Node(x=0, y=1, z=0)
    
    s1 = teefem.geom.Seg2(nodes = (n1,n2))
    s2 = teefem.geom.Seg2(nodes = (n2,n3))
    s3 = teefem.geom.Seg2(nodes = (n3,n1))
    
    mesh1 = teefem.geom.Mesh()
    mesh1.nodes = set([n1,n2,n3])
    mesh1.shapes = set([s1,s2,s3])
    mesh1.nset['GA'] = set([n1,n3])
    mesh1.nset['LOAD'] = set([n2])
    mesh1.elset['OM'] = set([s1,s2,s3])

    model = Bar2DModel(mesh = mesh1)
    mat = teefem.materials.elastic(E = 210.0e9, nu = 0.3)
    teefem.assign_material(elements = model.elset['OM'], material = mat)

    load = teefem.nodal_force(DY = -100.0e3)
    encastre = teefem.dirichlet_bc(encastre = True)
    
    teefem.assign_bc(nodes = model.nset['LOAD'], bc = load)
    teefem.assign_bc(nodes = model.nset['GA'], bc = encastre)

    carel = BarCharacteristic(area = 10e-4)

    teefem.assign_char(elements = model.elset['OM'], char = carel)

    model.static_solve()

    dymin = min([node.fields['DEPL']['DY'] for node in model.nodes])
    print("DY: %0.14E"%(dymin*1000))

if __name__ == '__main__':
    ex1()