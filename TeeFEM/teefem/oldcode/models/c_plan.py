# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 23:30:24 2012

@author: Jukka Aho


C_PLAN Plane stress model

"""

import teefem
from numpy import sqrt

class MEPLSE2(teefem.elements.Element1D):
    ''' Plane stress linear line element '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*2
        super(MEPLSE2, self).__init__(*args, **kwds)
        self.ipoints = [(0.0)]
        self.iweights = [(2.0)]

class MEPLSE3(teefem.elements.Element1D):
    ''' Plane stress quadratic line element '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*3
        super(MEPLSE3, self).__init__(*args, **kwds)
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf
        self.ipoints = [(-1/sqrt(3)),(1/sqrt(3))]
        self.iweights = [(1.0),(1.0)]

class MECPTR3(teefem.elements.Element2D):
    """ Plane stress linear triangle finite element CST
    
    EXAMPLES:     
    
    >>> n1 = Node(x=0.0, y=0.0)
    >>> n2 = Node(x=2.0, y=1.0)
    >>> n3 = Node(x=0.5, y=2.0)
    >>> tr3 = Tria3(name = 'tr3', nodes = (n1, n2, n3))
    >>> tr3.material = Elastic(E = 210e9, nu = 0.3)
    >>> el3 = MECPTR3(shape = tr3)
    >>> el3.plot(filename = 'img/MECPTR3.png')
    
    .. image:: img/MECPTR3.png

    >>> el3emat = el3.material_matrix
    >>> print el3emat
    [[  2.3e+11   6.9e+10   0.0e+00]
     [  6.9e+10   2.3e+11   0.0e+00]
     [  0.0e+00   0.0e+00   8.1e+10]]

    >>> el3bmat = el3.kinematic_matrix(1/3,1/3)
    >>> print el3bmat
    [[-0.3  0.   0.6  0.  -0.3  0. ]
     [ 0.  -0.4  0.  -0.1  0.   0.6]
     [-0.4 -0.3 -0.1  0.6  0.6 -0.3]]
     
    >>> el3kmat = el3.stiffness_matrix
    >>> print el3kmat
    [[  5.9e+10   3.2e+10  -5.7e+10  -3.0e+10  -1.6e+09  -2.5e+09]
     [  3.2e+10   8.6e+10  -2.4e+10   1.6e+09  -8.2e+09  -8.7e+10]
     [ -5.7e+10  -2.4e+10   1.3e+11  -2.1e+10  -7.7e+10   4.5e+10]
     [ -3.0e+10   1.6e+09  -2.1e+10   5.4e+10   5.1e+10  -5.6e+10]
     [ -1.6e+09  -8.2e+09  -7.7e+10   5.1e+10   7.9e+10  -4.3e+10]
     [ -2.5e+09  -8.7e+10   4.5e+10  -5.6e+10  -4.3e+10   1.4e+11]]
    """    
    
    def __init__(self, *args, **kwds):
        self.dimension = 2*3
        super(MECPTR3, self).__init__(*args, **kwds)
        self.has_stiffness = True
        self.ipoints = [(1/3, 1/3)]
        self.iweights = [(1/2)]
        
class MECPTR6(teefem.elements.Element2D):
    """ Plane stress quadratic triangular finite element 
    
    EXAMPLES:     
    
    >>> n1 = Node(x=0.0, y=0.0)
    >>> n2 = Node(x=2.0, y=1.0)
    >>> n3 = Node(x=0.5, y=2.0)
    >>> n4 = Node(x=1.0, y=0.0)
    >>> n5 = Node(x=1.0, y=1.0)
    >>> n6 = Node(x=0.2, y=1.0)
    >>> tr6 = Tria6(name = 'tr6', nodes = (n1, n2, n3, n4, n5, n6))
    >>> tr6.material = Elastic(E = 210e9, nu = 0.3)
    >>> el6 = MECPTR6(shape = tr6)
    >>> el6.plot(filename = 'img/MECPTR6.png')
    
    .. image:: img/MECPTR6.png

    >>> el6emat = el6.material_matrix
    >>> print el6emat
    [[  2.3e+11   6.9e+10   0.0e+00]
     [  6.9e+10   2.3e+11   0.0e+00]
     [  0.0e+00   0.0e+00   8.1e+10]]

    >>> el6bmat = el6.kinematic_matrix(1/3,1/3)
    >>> print el6bmat
    [[-0.2  0.   0.2  0.  -0.   0.   0.1  0.   0.7  0.  -0.8  0. ]
     [ 0.  -0.2  0.  -0.   0.   0.2  0.  -0.7  0.   0.6  0.   0.1]
     [-0.2 -0.2 -0.   0.2  0.2 -0.  -0.7  0.1  0.6  0.7  0.1 -0.8]]
     
    >>> el6kmat = el6.stiffness_matrix
    >>> print el6kmat
    [[  1.9e+11   7.2e+10   2.3e+10   1.3e+10   1.9e+10   1.7e+10  -1.5e+11  -2.0e+10  -9.5e+09  -2.5e+10  -7.4e+10  -5.7e+10]
     [  7.2e+10   1.4e+11   1.1e+10   8.4e+09   1.9e+10   2.7e+10  -1.2e+10  -4.4e+10  -2.5e+10  -1.2e+09  -6.4e+10  -1.3e+11]
     [  2.3e+10   1.1e+10   1.2e+11   4.6e+09   1.5e+10  -1.0e+10  -1.5e+10  -4.4e+10  -1.5e+11   4.8e+10   1.2e+10  -9.7e+09]
     [  1.3e+10   8.4e+09   4.6e+09   4.2e+10  -1.2e+10   7.1e+09  -5.1e+10  -2.1e+10   5.6e+10  -4.1e+10  -9.7e+09   3.7e+09]
     [  1.9e+10   1.9e+10   1.5e+10  -1.2e+10   3.6e+10   9.3e+09   9.9e+08  -2.2e+10   5.3e+09   7.0e+10  -7.5e+10  -6.5e+10]
     [  1.7e+10   2.7e+10  -1.0e+10   7.1e+09   9.3e+09   8.5e+10  -2.2e+10   1.1e+10   6.2e+10  -3.6e+10  -5.7e+10  -9.5e+10]
     [ -1.5e+11  -1.2e+10  -1.5e+10  -5.1e+10   9.9e+08  -2.2e+10   2.6e+11  -5.5e+10  -6.6e+10   4.1e+10  -2.8e+10   9.9e+10]
     [ -2.0e+10  -4.4e+10  -4.4e+10  -2.1e+10  -2.2e+10   1.1e+10  -5.5e+10   3.6e+11   4.1e+10  -2.7e+11   9.9e+10  -4.5e+10]
     [ -9.5e+09  -2.5e+10  -1.5e+11   5.6e+10   5.3e+09   6.2e+10  -6.6e+10   4.1e+10   5.9e+11  -9.6e+10  -3.7e+11  -3.8e+10]
     [ -2.5e+10  -1.2e+09   4.8e+10  -4.1e+10   7.0e+10  -3.6e+10   4.1e+10  -2.7e+11  -9.6e+10   4.4e+11  -3.8e+10  -9.4e+10]
     [ -7.4e+10  -6.4e+10   1.2e+10  -9.7e+09  -7.5e+10  -5.7e+10  -2.8e+10   9.9e+10  -3.7e+11  -3.8e+10   5.4e+11   7.1e+10]
     [ -5.7e+10  -1.3e+11  -9.7e+09   3.7e+09  -6.5e+10  -9.5e+10   9.9e+10  -4.5e+10  -3.8e+10  -9.4e+10   7.1e+10   3.6e+11]]
    """
    def __init__(self, *args, **kwds):
        self.dimension = 2*6
        super(MECPTR6, self).__init__(*args, **kwds)
        self.has_stiffness = True
        self.ipoints = [(1/6, 1/6),(2/3,1/6),(1/6,2/3)]
        self.iweights = [1/6, 1/6, 1/6]

class MECPQU4(teefem.elements.Element2D):
    """ Plane stress linear quadrangle finite element
    
    EXAMPLES:     
    
    >>> n1 = Node(x=1.0, y=0.0)
    >>> n2 = Node(x=2.5, y=0.5)
    >>> n3 = Node(x=2.0, y=2.0)
    >>> n4 = Node(x=-0.5, y=1.5)
    >>> qu4 = Quad4(name = 'MA1', nodes = (n1, n2, n3, n4))
    >>> qu4.material = Elastic(E = 100e9, nu = 0.5)
    >>> ma1 = MECPQU4(shape = qu4)
    >>> ma1.plot(filename = 'img/MECPQU4.png')
    
    .. image:: img/MECPQU4.png

    >>> ma1emat = ma1.material_matrix
    >>> print ma1emat
    [[  1.3e+11   6.7e+10   0.0e+00]
     [  6.7e+10   1.3e+11   0.0e+00]
     [  0.0e+00   0.0e+00   3.3e+10]]

    >>> ma1bmat = ma1.kinematic_matrix(0,0)
    >>> print ma1bmat
    [[-0.1  0.   0.3  0.   0.1  0.  -0.3  0. ]
     [ 0.  -0.4  0.  -0.1  0.   0.4  0.   0.1]
     [-0.4 -0.1 -0.1  0.3  0.4  0.1  0.1 -0.3]]
     
    >>> ma1kmat = ma1.stiffness_matrix
    >>> print ma1kmat
    [[  4.7e+10   2.3e+10  -3.0e+10  -1.1e+10  -1.9e+10  -2.0e+10   2.0e+09   8.6e+09]
     [  2.3e+10   1.1e+11  -2.8e+10  -2.7e+09  -2.0e+10  -7.0e+10   2.5e+10  -3.9e+10]
     [ -3.0e+10  -2.8e+10   6.1e+10  -1.2e+10  -1.5e+09   2.5e+10  -2.9e+10   1.5e+10]
     [ -1.1e+10  -2.7e+09  -1.2e+10   4.9e+10   8.3e+09  -4.4e+10   1.5e+10  -2.5e+09]
     [ -1.9e+10  -2.0e+10  -1.5e+09   8.3e+09   4.0e+10   2.2e+10  -1.9e+10  -1.0e+10]
     [ -2.0e+10  -7.0e+10   2.5e+10  -4.4e+10   2.2e+10   1.0e+11  -2.7e+10   1.3e+10]
     [  2.0e+09   2.5e+10  -2.9e+10   1.5e+10  -1.9e+10  -2.7e+10   4.7e+10  -1.4e+10]
     [  8.6e+09  -3.9e+10   1.5e+10  -2.5e+09  -1.0e+10   1.3e+10  -1.4e+10   2.8e+10]]
    
    """
    def __init__(self, *args, **kwds):
        self.dimension = 2*4
        super(MECPQU4, self).__init__(*args, **kwds)
        self.has_stiffness = True
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf s. 10
        a = 1/sqrt(3)
        self.ipoints = [(-a, -a),(a,-a),(a,a),(-a,a)] # Integroimispisteet
        self.iweights = [(1.0),(1.0),(1.0),(1.0)] # Painokertoimet ??????

class MECPQU8(teefem.elements.Element2D):
    ''' Plane stress quadratic 8-node Serendip finite element '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*8
        super(MECPQU8, self).__init__(*args, **kwds)
        self.has_stiffness = True
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf s. 10
        a = sqrt(3/5)
        self.ipoints = [(-a, -a),(a,-a),(a,a),(-a,a),(0,-a),(a,0),(0,a),(-a,0),(0,0)] # Integroimispisteet
        self.iweights = [(25/81),(25/81),(25/81),(25/81),(40/81),(40/81),(40/81),(40/81),(64/81)] # Painokertoimet

class MECPQU9(teefem.elements.Element2D):
    ''' Plane stress quadratic 9-node finite element '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*9
        super(MECPQU9, self).__init__(*args, **kwds)
        self.has_stiffness = True
        # http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf s. 10
        a = sqrt(3/5)
        self.ipoints = [(-a, -a),(a,-a),(a,a),(-a,a),(0,-a),(a,0),(0,a),(-a,0),(0,0)] # Integroimispisteet
        self.iweights = [(25/81),(25/81),(25/81),(25/81),(40/81),(40/81),(40/81),(40/81),(64/81)] # Painokertoimet

class C_PLAN(teefem.models.GenericModel):
    
    ''' Plane stress modelisation '''
    
    def __init__(self, *args, **kwds):
        super(C_PLAN, self).__init__(*args, **kwds)
        mapping = {
            'Seg2': MEPLSE2,
            'Seg3': MEPLSE3,
            'Tria3': MECPTR3,
            'Tria6': MECPTR6,
            'Quad4': MECPQU4,
            'Quad8': MECPQU8,
            'Quad9': MECPQU9,
        }
        # Create elements according element mapping
        for (name, net) in self.mesh.nets.iteritems():
            try:
                shape = net.__class__.__name__
                self.elements[name] = mapping[shape](shape = net)
            except KeyError:
                logging.debug("No corresponding finite element to net {0} in modelisation C_PLAN".format(shape))
        self.group_ma = {}
        for (k,v) in self.mesh.group_ma.iteritems():
            try: 
                self.group_ma[k] = set([self.elements[net.name] for net in v])
                logging.debug("New finite element group: {0}".format(k))
            except KeyError:
                logging.debug("Unable to create finite element group for group_ma {0}".format(k))
    
    def assembly_stiffness_matrix(self):
        ''' Return C_PLAN stiffness matrix '''
        # 0) Size of stiffness matrix
        ijvsize = 0
        kdim = 0
        for element in self.elements.itervalues():
            if not element.has_stiffness:
                continue
            ijvsize += element.dimension**2
            kdim += element.dimension
        # 1) NUME_DDL (assign position to global stiffness matrix)
        idx = 0
        for node in self.nodes.itervalues():
            node.gdof = xrange(idx*2,(idx+1)*2)
            if hasattr(node, 'bc'):
                if node.bc.has_key('dx'):
                    ijvsize += 2
                if node.bc.has_key('dy'):
                    ijvsize += 2
            idx += 1
        ijv = IJV(ijvsize)
        # 2) Assembly element local stiffness matrix to global stiffness matrix
        for element in self.elements.itervalues():
            if not element.has_stiffness:
                continue
            k_loc = element.stiffness_matrix
            gdof = array([node.gdof for node in element.nodes]).flatten()
            for I in xrange(element.dimension):
                for J in xrange(element.dimension):
                    ijv.add(gdof[I], gdof[J], k_loc[I,J])
        # 3) Assign boundary conditions to stiffness matrix with Lagrange multipliers
        alp = 1e12
        lidx = len(self.nodes)*2
        for node in self.nodes.itervalues():
            if hasattr(node, 'bc'):
                if node.bc.has_key('dx'):
                    ijv.add(node.gdof[0],lidx,alp)
                    ijv.add(lidx,node.gdof[0],alp)
                    lidx += 1
                if node.bc.has_key('dy'):
                    ijv.add(node.gdof[1],lidx,alp)
                    ijv.add(lidx,node.gdof[1],alp)
                    lidx += 1
        self.total_dimension = lidx # Total dimension of stiffness matrix
        return ijv.tocoo()
    
    def assembly_force_vector(self):
        ''' Return C_PLAN force vector '''
        R = np.zeros(self.total_dimension)
        for element in self.elements.itervalues():
            if not hasattr(element, 'bc'):
                continue
            bc = element.bc
            element.qn = lambda k: bc['pressure']
            element.qt = lambda k: bc['shear']
            fq = element.get_force_vector()
            gdof = array([node.gdof for node in element.nodes]).flatten()
            for I in xrange(element.dimension):
                R[gdof[I]] += fq[I]
        lidx = len(self.nodes)*2
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
        
    def static_solve(self):
        logging.debug("Creating stiffness matrix")
        K = self.get_stiffness_matrix().tocsr()
        logging.debug("Creating force vector")
        R = self.get_force_vector()
        logging.debug("Solving")
        U = spsolve(K,R)
        logging.debug("Updating fields")
        for node in self.nodes.itervalues():
            node.update_field(
                field = 'DEPL', 
                params = ('DX','DY'), 
                values = (U[node.gdof[0]], U[node.gdof[1]]),
                )
