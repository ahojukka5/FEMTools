# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 12:04:59 2012

@author: Jukka

"""

import teefem
from teefem import log, cache
#import numpy as np
from numpy import array, zeros

def assign_material(**kwds):
    elements = kwds.get('elements')
    material = kwds.get('material')
    log.info('Assigning material to {0} elements...'.format(len(elements)))
    for element in elements:
        #element.material = material
        try:
            element.assign_material(material)
        except AttributeError as err:
            log.error(err)
    return True

def assign_bc(**kwds):
    bc = kwds.get('bc')
    for element in kwds.get('elements', []):
        try:
            element.assign_boundary_condition(bc)
        except AttributeError as err:
            log.error(err)
    for node in kwds.get('nodes', []):
        node.assign_boundary_condition(bc)
    return True

def assign_char(**kwds):
    char = kwds.get('char')
    for element in kwds.get('elements', []):
        try:
            element.assign_char(char)
        except AttributeError as err:
            log.error(err)
    return True

class Model(object):

    ''' General mechanical model '''

    def __init__(self, *args, **kwds):
        self.mesh = kwds['mesh']
        self.alp = 1e12
        teefem.log.debug("Model created")

    @teefem.deprecated
    def assign_bc(self,*args,**kwds):
        ''' Assigns boundary conditions to elements / nodes '''
        if kwds.has_key('elem_group'):
            for element in self.elset[kwds['elem_group']]:
                element.bc = kwds
        if kwds.has_key('node_group'):
            try:
                for node in self.nset[kwds['node_group']]:
                    node.bc = kwds
            except KeyError:
                log("Error: node group {0} does not exist!".format(kwds['node_group']))
    
    @property
    @cache
    def stiff_elements(self):
        return set([element for element in self.elements if hasattr(element, 'stiffness_matrix')])

    @property
    @cache
    def load_elements(self):
        return set([element for element in self.elements if hasattr(element, 'force_vector')])

    def init(self):
        mapping = self.mapping

        # Create elements according element mapping
        # Add element to elements set, add nodes to nodes set
        self.elements = set()
        self.nodes = set()
        for (name, shape) in self.mesh.nets.iteritems():
            try:
                shapename = shape.__class__.__name__
                element = mapping[shapename](geom = shape)
                self.elements.add(element)
                for node in element.nodes:
                    self.nodes.add(node)
            except KeyError:
                log.warning("No corresponding finite element to shape {0}".format(shape))
        
        # Create element groups and node groups
        self.create_node_groups()
        self.create_element_groups()
    
    def create_node_groups(self):
        log.info('Creating node groups for model')
        self.nset = {}
        for node in self.nodes:
            for ngroup in node.groups:
                try:
                    self.nset[ngroup].add(node)
                except KeyError:
                    log.info('New node group: {0}'.format(ngroup))
                    self.nset[ngroup] = set([node])
    
    def create_element_groups(self):
        log.info('Creating element groups for model')
        self.elset = {}
        for element in self.elements:
            for egroup in element.groups:
                try:
                    self.elset[egroup].add(element)
                except KeyError:
                    log.info('New element group: {0}'.format(egroup))
                    self.elset[egroup] = set([element])

#        for (k,v) in self.mesh.group_ma.iteritems():
#            try: 
#                self.group_ma[k] = set([self.elements[net.name] for net in v])
#                log.info("New finite element group: {0}".format(k))
#            except KeyError:
#                log.warning("Unable to create finite element group for group_ma {0}".format(k))
    
    def _init_stiffness_matrix(self):
        ''' Initialize stiffness matrix '''
        
        nodedim = self.nodedim
        self.alp = 1e12
        self.lidx = len(self.nodes)*nodedim

        for node in self.nodes:
            node.dim = nodedim

        idx = 0

        for node in self.nodes:
            node.update_boundary_conditions()
            node.gdof = xrange(idx,idx+node.dim)
            idx += node.dim

        self.ijv = teefem.common.IJV()

    def _populate_ijv(self):
        ''' Populate IJV coordinate system '''
        for element in self.stiff_elements:
            try:
                k_loc = element.stiffness_matrix
                gdof = array([node.gdof for node in element.nodes]).flatten()
                for I in xrange(element.dimension):
                    for J in xrange(element.dimension):
                        self.ijv.add(gdof[I], gdof[J], k_loc[I,J])
            except AttributeError as err:
                log.error(err)

    def _assign_boundary_conditions_to_stiffness_matrix(self):
        ''' Assign boundary conditions to IJV matrix '''
        alp = self.alp
        ijv = self.ijv
        lidx = len(self.nodes)*self.nodedim
        for node in self.nodes:
            for j in range(self.nodedim):
                dof = self.nodedofs[j]
                if node.hasbc(dof):
                    ijv.add(node.gdof[j],lidx,alp)
                    ijv.add(lidx,node.gdof[j],alp)
                    lidx += 1
        self.total_dimension = lidx # Total dimension of stiffness matrix

#    def _assign_boundary_conditions_to_stiffness_matrix(self):
#        ''' Assign boundary conditions to IJV coordinate system '''
#        alp = self.alp
#        ijv = self.ijv
#        lidx = len(self.nodes)*self.nodedim
#        for node in self.nodes:
#            if hasattr(node, 'bc'):
#                if node.bc.has_key('dx'):
#                    ijv.add(node.gdof[0],lidx,alp)
#                    ijv.add(lidx,node.gdof[0],alp)
#                    lidx += 1
#                if node.bc.has_key('dy'):
#                    ijv.add(node.gdof[1],lidx,alp)
#                    ijv.add(lidx,node.gdof[1],alp)
#                    lidx += 1
#        self.total_dimension = lidx # Total dimension of stiffness matrix

    def assemble_stiffness_matrix(self):
        ''' Assemble stiffness matrix '''
        self._init_stiffness_matrix()
        self._populate_ijv()
        self._assign_boundary_conditions_to_stiffness_matrix()
        return self.ijv.tocoo()

    def _init_force_vector(self):    
        ''' Initialize force vector '''
        self.R = zeros(self.total_dimension)

    def _populate_R(self):    
        R = self.R
        for element in self.load_elements:
#            if not hasattr(element, 'bc'):
#                continue
#            bc = element.bc
#            element.qn = lambda k: bc['pressure']
#            element.qt = lambda k: bc['shear']
            fq = element.force_vector
            gdof = array([node.gdof for node in element.nodes]).flatten()
            for I in xrange(element.dimension):
                R[gdof[I]] += fq[I]

#    def _assign_boundary_conditions_to_force_vector(self):
#        lidx = len(self.nodes)*self.nodedim
#        R = self.R
#        alp = self.alp
#        for node in self.nodes:
#            if hasattr(node, 'bc'):
#                if node.bc.has_key('fx'):
#                    R[node.gdof[0]] += node.bc['fx']
#                if node.bc.has_key('fy'):
#                    R[node.gdof[1]] += node.bc['fy']
#                if node.bc.has_key('dx'):
#                    R[lidx] = node.bc['dx']*alp
#                    lidx += 1
#                if node.bc.has_key('dy'):
#                    R[lidx] = node.bc['dy']*alp
#                    lidx += 1
#        return R

    def _assign_boundary_conditions_to_force_vector(self):
        lidx = len(self.nodes)*self.nodedim
        R = self.R
        alp = self.alp
        for node in self.nodes:
            for j in range(self.nodedim):
                dof = self.nodedofs[j]
                if node.hasbc(dof):
                    R[lidx] = node.valbc(dof)*alp
                    lidx += 1
                load = self.nodeloads[j]
                if node.hasload(load):
                    R[node.gdof[j]] += node.valload(load)
        return R
        
    def assemble_force_vector(self):
        ''' Assemble force vector '''
        self._init_force_vector()
        self._populate_R()
        self._assign_boundary_conditions_to_force_vector()
        return self.R

    def static_solve(self, **kwds):
        log.info("Creating stiffness matrix")
        K = self.assemble_stiffness_matrix().tocsr()
        log.info("Creating force vector")
        R = self.assemble_force_vector()
        log.info("Solving")
        U = teefem.common.solve(K,R)
        log.info("Updating fields")
        # log.debug(U)
        for element in self.elements:
            element.update(U)
