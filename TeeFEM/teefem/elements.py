# -*- coding: utf-8 -*-
"""
Created on Thu Mar 08 00:03:31 2012

@author: Jukka Aho

#####################################
### Elements.py - FINITE ELEMENTS ###
#####################################
"""

from __future__ import division
# import sys
# import numpy as np
from numpy import matrix, zeros
from teefem import cache

class Element(object):

    """ General finite element class """    
    
    def __init__(self, *args, **kwds):
        self.boundary_conditions = set()
        # Perusvaatimus, jokaisella elementillä on:
        # 1) Geometria (geom,shape)
        # 2) Solmupisteet (geometrian perusteella)
        self.geom = kwds['geom']
        self.nodes = self.geom.nodes # Isoparametric
        self.groups = self.geom.groups
#        self.has_stiffness = False
#        self.material = self.geom.material
        
    @property
    @cache
    def detJ(self):
        return self.geom.detJ

    @property
    def J(self):
        return self.geom.J
    
    @property
    def invJ(self):
        return self.geom.invJ
    
    @property
    def area(self):
        return self.geom.area

    def update(self, U):
        u = zeros(self.dimension)
        nodedim = len(self.degrees_of_freedom)
        # Tämä kierrätetään ihan turhaan elementin kautta
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

    def update_modal(valsy, vecsy):
        pass

    @property
    def status(self):
        val = ''
        val += 'Element type: {0}\n'.format(self.__class__.__name__)
        val += 'Element geometry type: {0}\n'.format(self.geom.__class__.__name__)
        val += 'No of DOFs: {0}\n\n'.format(self.dimension)
        val += 'Geometry:\n'
        val += 'Number of nodes: {0}\n'.format(len(self.geom.nodes))
        val += 'Node coordinates:\n'
        val += '  #        x        y        z\n'
        for i in range(len(self.geom.nodes)):
            val += '{0:3d} {node.x:7.2E} {node.x:7.2E} {node.x:7.2E}\n'.format(i, node = self.geom.nodes[i])
        return val

class Element1D(Element):
    
    ''' General functions for 1D finite elements '''
    
    def __init__(self, *args, **kwds):
        super(Element1D, self).__init__(*args, **kwds)

#    def get_force_vector(self):
#        shape = self.geom
#        Rq = zeros(self.dimension)
#        qx = array([W * self.qn(x) * shape.N(x,0) * shape.J12(x,0) for (W,x) in zip(self.iweights,self.ipoints)])
#        qy = array([W * self.qn(x) * shape.N(x,0) * shape.J11(x,0) for (W,x) in zip(self.iweights,self.ipoints)])
#        Rq[1::2] += np.sum(qy,axis=0)
#        Rq[0::2] += np.sum(qx,axis=0)
#        return Rq

class Element2D(Element):
    
    ''' General 2D finite element class '''

    def __init__(self, *args, **kwds):
        super(Element2D, self).__init__(*args, **kwds)


#    def kinematic_matrix(self,*ke):
#        r""" Returns plane kinematic matrix B 
#
#        .. math ::
#        
#            \begin{align}
#            \mathbf{B} & =\begin{bmatrix}\mathbf{B}_{1} & \mathbf{B}_{2} & \mathbf{B}_{3} & \mathbf{B}_{4} & \mathbf{B}_{5} & \mathbf{B}_{6}\end{bmatrix} & \mathbf{B_{i}} & =\begin{bmatrix}\frac{\partial N_{i}}{\partial x} & 0\\
#            0 & \frac{\partial N_{i}}{\partial y}\\
#            \frac{\partial N_{i}}{\partial y} & \frac{\partial N_{i}}{\partial x}
#            \end{bmatrix}
#            \end{align}
#        
#        INPUT:
#        - (``k``,``e``)
#        
#        OUTPUT:
#        - 3x2m matrix
#        
#        """
#        dNdk = self.dNdk(*ke)
#        dNde = self.dNde(*ke)
#        dNxy = self.invJ(*ke) * matrix([dNdk,dNde])
#        dNdx = dNxy[0]
#        dNdy = dNxy[1]
#        B = matrix(zeros((3,self.dimension)))
#        B[0,0::2] = dNdx
#        B[1,1::2] = dNdy
#        B[2,0::2] = dNdy
#        B[2,1::2] = dNdx
#        return B
            
#    @property
#    def material_matrix(self):
#        r""" Returns plane material matrix E
#        
#        .. math ::
#            
#            \mathbf{E}=\frac{E}{1-\nu^{2}}\begin{bmatrix}1 & \nu & 0\\
#            \nu & 1 & 0\\
#            0 & 0 & \frac{1-\nu}{2}
#            \end{bmatrix}
#            
#        """
#        try:
#            return self._E
#        except AttributeError:
#            EE = self.geom.material.E
#            nu = self.geom.material.nu
#            self._E = EE/(1-nu**2)*matrix([[1,nu,0],[nu,1,0],[0,0,(1-nu)/2]])
#            return self._E

#    @property
#    def stiffness_matrix(self):
#        geom = self.geom
#        detJ = geom.detJ
#        B = self.kinematic_matrix
#        E = self.material_matrix
#        b = array([W * B(*ke).T * E * B(*ke) * detJ(*ke) for (W,ke) in zip(self.iweights, self.ipoints)])
#        return np.sum(b,axis=0)
#
    def plot(self, *args, **kwds):
        ''' Plots geometry and Gauss points'''
        self.geom.plot()
        x = [self.geom.x(*p) for p in self.ipoints]
        y = [self.geom.y(*p) for p in self.ipoints]
        plt.scatter(x,y,marker='o',c='r',s=10)
        for (xi,yi,Wi) in zip(x,y,self.iweights): 
            plt.annotate('%0.2f'%(Wi+1), xy=(xi,yi), xytext=(0,10), 
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='red', alpha=0.2))
        plt.title(self.__class__.__name__)
        if kwds.has_key('filename'):
            plt.savefig(**kwds)


