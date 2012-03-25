# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 10:38:43 2012

@author: Jukka
"""


class BoundaryCondition(object):
    pass

class PressureBoundaryCondition(BoundaryCondition):
    '''
    INPUT KEYWORDS
    
        pressure : function
        f.e. pressure = lambda k,e: -100e3
    '''
    def __init__(self, **kwds):
        self.pressure = kwds.get('pressure')

class NodalForce(BoundaryCondition):

    def __init__(self, **kwds):
        self._bc = {}
        if kwds.has_key('fx'):
            self._bc['fx'] = kwds['fx']
        if kwds.has_key('fy'):
            self._bc['fy'] = kwds['fy']
        if kwds.has_key('fz'):
            self._bc['fz'] = kwds['fz']
        if kwds.has_key('mx'):
            self._bc['mx'] = kwds['mx']
        if kwds.has_key('my'):
            self._bc['my'] = kwds['my']
        if kwds.has_key('mz'):
            self._bc['mz'] = kwds['mz']
    def exist(self,dof):
        return self._bc.has_key(dof)
    def val(self,dof):
        return self._bc[dof]

class DirichletBoundaryCondition(BoundaryCondition):
    def __init__(self, **kwds):
        self._bc = {}
        if kwds.has_key('dx'):
            self._bc['dx'] = kwds['dx']
        if kwds.has_key('dy'):
            self._bc['dy'] = kwds['dy']
        if kwds.has_key('dz'):
            self._bc['dz'] = kwds['dz']
        if kwds.has_key('drx'):
            self._bc['drx'] = kwds['drx']
        if kwds.has_key('dry'):
            self._bc['dry'] = kwds['dry']
        if kwds.has_key('drz'):
            self._bc['drz'] = kwds['drz']
        if kwds.get('encastre',False):
            self._bc['dx'] = 0
            self._bc['dy'] = 0
            self._bc['dz'] = 0
            self._bc['drx'] = 0
            self._bc['dry'] = 0
            self._bc['drz'] = 0
            self._bc['gx'] = 0
    def exist(self,dof):
        return self._bc.has_key(dof)
    def val(self,dof):
        return self._bc[dof]

pressure_bc = PressureBoundaryCondition
dirichlet_bc = DirichletBoundaryCondition
nodal_force = NodalForce

def test1():
    bc = DirichletBoundaryCondition(dx = -1, dy = 0, drz = 1, foo = 0)
    print bc._bc
    print bc.exist('dry')
    print bc.exist('drz')
    print bc.val('drz')
    print bc.val('dx')
    bc2 = DirichletBoundaryCondition(encastre = True)
    print bc2._bc

if __name__ == '__main__':
    test1()