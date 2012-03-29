# -*- coding: utf-8 -*-
"""
Created on Sun Mar 11 12:04:59 2012

@author: Jukka Aho

Yleisiä mekaaniseen malliin kuuluvia ominaisuuksia

27.3.2012. DRY.

"""

import teefem
from teefem import log, cache
#import numpy as np
from numpy import array, zeros
import numpy as np
import scipy

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
        self.penalty_term = 1e15
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
        return set([element for element in self.elements if element.has_stiffness_matrix])

    @property
    @cache
    def geometric_stiff_elements(self):
        return set([element for element in self.elements if element.has_geometric_stiffness_matrix])

    @property
    @cache
    def load_elements(self):
        return set([element for element in self.elements if element.has_load_vector])

    @property
    @cache
    def mass_elements(self):
        return set([element for element in self.elements if element.has_mass_matrix])

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
        self._nset = {}
        for node in self.nodes:
            for ngroup in node.groups:
                try:
                    self._nset[ngroup].add(node)
                except KeyError:
                    log.info('New node group: {0}'.format(ngroup))
                    self._nset[ngroup] = set([node])
    
    def create_element_groups(self):
        log.info('Creating element groups for model')
        self._elset = {}
        for element in self.elements:
            for egroup in element.groups:
                try:
                    self._elset[egroup].add(element)
                except KeyError:
                    log.info('New element group: {0}'.format(egroup))
                    self._elset[egroup] = set([element])

    @cache
    def nset(self,*args):
        nset = set()
        for seti in args:
            nodes = self._nset[seti]
            for node in nodes:
                nset.add(node)
        return nset

    @cache
    def elset(self,*args):
        elset = set()
        for sete in args:
            elements = self._elset[sete]
            for element in elements:
                elset.add(element)
        return elset

#        for (k,v) in self.mesh.group_ma.iteritems():
#            try: 
#                self.group_ma[k] = set([self.elements[net.name] for net in v])
#                log.info("New finite element group: {0}".format(k))
#            except KeyError:
#                log.warning("Unable to create finite element group for group_ma {0}".format(k))

###############################################################################

    
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
        return teefem.common.IJV()

    def _populate_stiffness_matrix(self,ijv):
        ''' Populate IJV coordinate system '''
        for element in self.stiff_elements:
            try:
                k_loc = element.stiffness_matrix
                gdof = array([node.gdof for node in element.nodes]).flatten()
                for I in xrange(element.dimension):
                    for J in xrange(element.dimension):
                        ijv.add(gdof[I], gdof[J], k_loc[I,J])
            except AttributeError as err:
                log.error(err)
        return ijv

    def _assign_lagrange_boundary_conditions_to_stiffness_matrix(self,ijv):
        ''' Assign boundary conditions to IJV matrix, Lagrange method '''
        alp = self.alp
        lidx = len(self.nodes)*self.nodedim
        for node in self.nodes:
            for j in range(self.nodedim):
                dof = self.nodedofs[j]
                if node.hasbc(dof):
                    ijv.add(node.gdof[j],lidx,alp)
                    ijv.add(lidx,node.gdof[j],alp)
                    lidx += 1
        self.total_dimension = lidx # Total dimension of stiffness matrix
        return ijv

    def _assign_penalty_boundary_conditions_to_stiffness_matrix(self):
        ''' Assign boundary conditions to IJV matrix, penalty method '''
        pterm = self.penalty_term
        ijv = self.ijv
        for node in self.nodes:
            for j in range(self.nodedim):
                dof = self.nodedofs[j]
                if node.hasbc(dof):
                    ijv.add(node.gdof[j],node.gdof[j],pterm)
        # Total dimension of stiffness matrix
        self.total_dimension = len(self.nodes)*self.nodedim 


    def assemble_stiffness_matrix(self, **kwds):
        ''' Assemble stiffness matrix '''
        method = kwds.get('method','lagrange')
        ijv = self._init_stiffness_matrix()
        ijv = self._populate_stiffness_matrix(ijv)
        if method == 'lagrange':
            ijv = self._assign_lagrange_boundary_conditions_to_stiffness_matrix(ijv)
        if method == 'penalty':
            ijv = self._assign_penalty_boundary_conditions_to_stiffness_matrix(ijv)
        return ijv.tocoo()

###############################################################################

    def _init_mass_matrix(self):
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

    def _populate_mass_matrix(self):
        ''' Populate IJV coordinate system '''
        for element in self.mass_elements:
            try:
                k_loc = element.mass_matrix
                gdof = array([node.gdof for node in element.nodes]).flatten()
                for I in xrange(element.dimension):
                    for J in xrange(element.dimension):
                        self.ijv.add(gdof[I], gdof[J], k_loc[I,J])
            except AttributeError as err:
                log.error(err)

#
#    def _assign_penalty_boundary_conditions_to_mass_matrix(self):
#        ''' Assign boundary conditions to IJV matrix '''
#        alp = self.alp
#        ijv = self.ijv
#        lidx = len(self.nodes)*self.nodedim
#        for node in self.nodes:
#            for j in range(self.nodedim):
#                dof = self.nodedofs[j]
#                if node.hasbc(dof):
#                    ijv.add(node.gdof[j],lidx,alp)
#                    ijv.add(lidx,node.gdof[j],alp)
#                    lidx += 1
#        self.total_dimension = lidx # Total dimension of stiffness matrix

    def assemble_mass_matrix(self, **kwds):
        ''' Assemble mass matrix '''
        self._init_mass_matrix()
        self._populate_mass_matrix()
#        self._assign_boundary_conditions_to_mass_matrix()
        return self.ijv.tocoo()


###############################################################################

    def _init_geometric_stiffness_matrix(self):
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
        return teefem.common.IJV()

    def _populate_geometric_stiffness_matrix(self, ijv):
        ''' Populate IJV coordinate system '''
        for element in self.geometric_stiff_elements:
            try:
                k_loc = element.geometric_stiffness_matrix
                gdof = array([node.gdof for node in element.nodes]).flatten()
                for I in xrange(element.dimension):
                    for J in xrange(element.dimension):
                        ijv.add(gdof[I], gdof[J], k_loc[I,J])
            except AttributeError as err:
                log.error(err)
        return ijv

    def assemble_geometric_stiffness_matrix(self, **kwds):
        ''' Assemble geometric stiffness matrix Kg'''
        ijv = self._init_geometric_stiffness_matrix()
        ijv = self._populate_geometric_stiffness_matrix(ijv)
#        self._assign_boundary_conditions_to_mass_matrix()
        return ijv.tocoo()

###############################################################################        

    def _init_force_vector(self):    
        ''' Initialize force vector '''
        return zeros(self.total_dimension)

    def _populate_R(self, R):    
        for element in self.load_elements:
            fq = element.force_vector
            gdof = array([node.gdof for node in element.nodes]).flatten()
            for I in xrange(element.dimension):
                R[gdof[I]] += fq[I]
        return R

    def _assign_lagrange_boundary_conditions_to_force_vector(self,R):
        lidx = len(self.nodes)*self.nodedim
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
        
    def assemble_force_vector(self, **kwds):
        ''' Assemble force vector '''
        method = kwds.get('method','lagrange')
        R = self._init_force_vector()
        R = self._populate_R(R)
        if method == 'lagrange':
            R = self._assign_lagrange_boundary_conditions_to_force_vector(R)
        if method == 'penalty':
            R = self._assign_penalty_boundary_conditions_to_force_vector(R)
        return R

###############################################################################

    def static_solve(self, **kwds):
        log.info("Creating stiffness matrix")
        K = self.assemble_stiffness_matrix(method = 'lagrange').tocsr()
        log.info("Creating force vector")
        R = self.assemble_force_vector(method = 'lagrange')
        log.info("Solving")
        U = teefem.common.solve(K,R)
        if kwds.get('export_matrices',False):
            np.savetxt('K.txt', K.todense(), fmt='%+05.2E')
            np.savetxt('R.txt', R, fmt='%+05.2E')
            np.savetxt('U.txt', U, fmt='%+05.2E')        
        log.info("Updating fields")
        for element in self.elements:
            element.update(U)

###############################################################################

    def modal_solve(self, **kwds):
        
        ''' Modal analysis '''
        
        n = kwds.get('n_modes', 1)
        log.info("Creating stiffness matrix")
        K = self.assemble_stiffness_matrix(method = 'none').todense()
        log.info("Creating mass matrix")
        M = self.assemble_mass_matrix(method = 'none').todense()
        
        # Eliminate dofs w boundary conditions
        
        from teefem import remdof
        
        doflist = []
        for node in self.nodes:
            for j in range(self.nodedim):
                dof = self.nodedofs[j]
                if node.hasbc(dof):
                    doflist.append(node.gdof[j])
        
        K = remdof.removedofs(K,doflist)
        M = remdof.removedofs(M,doflist)
        
        if kwds.get('export_matrices',False):
            np.savetxt('K.txt', K, fmt='%+05.2E')
            np.savetxt('M.txt', M, fmt='%+05.2E')

        from scipy.linalg import eig
        valsy,vecsy = eig(K,b=M)
        
        print("w (Y): {0}".format(np.sort(np.sqrt(valsy.real))))
        
        #print("freq (Y): {0}".format(np.sort(np.sqrt(valsy.real)/(2*np.pi))))
        if kwds.get('export_matrices',False):
            np.savetxt('vals1.txt', np.sqrt(valsy.real), fmt='%+05.2E')
            np.savetxt('vecs1.txt', vecsy, fmt='%+05.2E')

#        log.info("Solving")
#        w,v = scipy.sparse.linalg.eigs(K, k=n, M=M)
#        print w
##
#        if kwds.get('export_matrices',False):
#            np.savetxt('vals2.txt', w, fmt='%+05.2E')
#            np.savetxt('vecs2.txt', v, fmt='%+05.2E')

        # Tämä on nyt varsin huonosti tehty. Kentät pitäisi päivittää jotenkin
        # järkevämpään tietorakenteeseen.

        log.info("Updating primary field")
        
        for node in self.nodes:
            node.modalvals = []
            node.modalvecs = []

        modal_indices = np.argsort(np.sqrt(valsy.real))
        for modal_idx in modal_indices[0:n]:
            val = np.sqrt(valsy[modal_idx].real)
            print val
            vec = vecsy[:,modal_idx]
            newvec = remdof.insertdofs(vec,doflist)
            for node in self.nodes:
                dof = node.gdof[0]
                node.modalvals.append(val)
                node.modalvecs.append(newvec[dof])

###############################################################################

    def buckling_solve(self, **kwds):
        
        ''' Buckling analysis '''
        
        n = kwds.get('n_modes', 1)
        log.info("Creating stiffness matrix")
        K = self.assemble_stiffness_matrix(method = 'none').todense()
        log.info("Creating geometric stiffness matrix")
        S = self.assemble_geometric_stiffness_matrix(method = 'none').todense()
        
        # Eliminate dofs w boundary conditions
        
        from teefem import remdof
        
        doflist = []
        for node in self.nodes:
            for j in range(self.nodedim):
                dof = self.nodedofs[j]
                if node.hasbc(dof):
                    doflist.append(node.gdof[j])
        
        K = remdof.removedofs(K,doflist)
        S = remdof.removedofs(S,doflist)
        
        if kwds.get('export_matrices',False):
            np.savetxt('K.txt', K, fmt='%+05.2E')
            np.savetxt('S.txt', S, fmt='%+05.2E')

        from scipy.linalg import eig
        valsy,vecsy = eig(K,b=S)
        
        print("Lambda : {0}".format(np.sort(np.sqrt(valsy.real))))
        
        #print("freq (Y): {0}".format(np.sort(np.sqrt(valsy.real)/(2*np.pi))))
        if kwds.get('export_matrices',False):
            np.savetxt('vals1.txt', np.sqrt(valsy.real), fmt='%+05.2E')
            np.savetxt('vecs1.txt', vecsy, fmt='%+05.2E')

#        log.info("Solving")
#        w,v = scipy.sparse.linalg.eigs(K, k=n, M=M)
#        print w
##
#        if kwds.get('export_matrices',False):
#            np.savetxt('vals2.txt', w, fmt='%+05.2E')
#            np.savetxt('vecs2.txt', v, fmt='%+05.2E')

        log.info("Updating primary field")
        
        # Tämä on nyt varsin huonosti tehty. Kentät pitäisi päivittää jotenkin
        # järkevämpään tietorakenteeseen
        
        for node in self.nodes:
            node.bucklingvals = []
            node.bucklingvecs = []

        buckling_indices = np.argsort(np.sqrt(valsy.real))
        for buckling_idx in buckling_indices[0:n]:
            val = np.sqrt(valsy[buckling_idx].real)
            print val
            vec = vecsy[:,buckling_idx]
            newvec = remdof.insertdofs(vec,doflist)
            for node in self.nodes:
                dof = node.gdof[0]
                node.bucklingvals.append(val)
                node.bucklingvecs.append(newvec[dof])