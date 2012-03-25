# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:13:31 2012

@author: Jukka Aho
"""

########################
### nets.py - SHAPES ###
########################

from numpy import matrix, concatenate, linspace, array
from numpy.linalg import det, inv
import matplotlib.pylab as plt

class Node(object):
    """ Node object
    
    INPUT:
        - ``x``
        - ``y``
        - ``z``
        
    OUTPUT:
        - Node object

    EXAMPLES:
    
    >>> n = Node(x = 1, y = 2, z = 3)
    >>> print n
    Node : (1.00,2.00,3.00)

    """
    def __init__(self, x=0, y=0, z=0, **kwds):
        self.x = x
        self.y = y
        self.z = z
        self.fields = {}
        
    def __repr__(self):
        return('{0} : ({1.x:0.2f},{1.y:0.2f},{1.z:0.2f})'.format(self.__class__.__name__, self))
        
    def update_field(self, *args, **kwds):
        field = kwds['field'].upper()
        if not self.fields.has_key(field):
            self.fields[field] = {}
        for (p,v) in zip(kwds['params'], kwds['values']):
            self.fields[kwds['field']][p] = v


class Net(object):
    def __init__(self, *args, **kwds):
        ''' General Net object
        
        EXAMPLES:        
        
        >>> node1 = Node(x=1,y=2,z=3)
        >>> node2 = Node(x=2,y=3,z=4)
        >>> node3 = Node(x=3,y=4,z=5)
        >>> net1 = Net(name = 'net1', nodes = [node1, node2, node3])
        
        '''
        self.name = kwds['name']
        self.nodes = kwds['nodes']
        self.xcoords = [n.x for n in self.nodes]
        self.ycoords = [n.y for n in self.nodes]
        self.x = lambda *ke: sum(N*nodes.x for (N,nodes) in zip(self.N(*ke),self.nodes))
        self.y = lambda *ke: sum(N*nodes.y for (N,nodes) in zip(self.N(*ke),self.nodes))

class Net1D(Net):
    def __init__(self, *args, **kwds):
        super(Net1D, self).__init__(*args, **kwds) 
        self.J11 = lambda *ke: sum(dNdk*nodes.x for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))
        self.J12 = lambda *ke: sum(dNdk*nodes.y for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))

class Net2D(Net):
    
    """ General 2d shapeset """
    
    def __init__(self, *args, **kwds):
        super(Net2D, self).__init__(*args, **kwds)
        self._J = {}
        self._detJ = {}
        self._invJ = {}

    def J(self,*ke):
        """ Returns 2d Jacobian matrix 
        
        EXAMPLES:        
        """
        try:
            return self._J[ke]
        except KeyError:
            J11 = sum(dNdk*nodes.x for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))
            J12 = sum(dNdk*nodes.y for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))
            J21 = sum(dNde*nodes.x for (dNde,nodes) in zip(self.dNde(*ke),self.nodes))
            J22 = sum(dNde*nodes.y for (dNde,nodes) in zip(self.dNde(*ke),self.nodes))
            self._J[ke] = matrix([[J11,J12],[J21,J22]])
            return self._J[ke]
    
    def invJ(self,*ke):
        ''' Returns inverse Jacobian matrix '''
        try:
            return self._invJ[ke]
        except KeyError:
            self._invJ[ke] = inv(self.J(*ke))
            return self._invJ[ke]
    
    def detJ(self,*ke):
        ''' Returns Jacobian '''
        try:
            return self._detJ[ke]
        except KeyError:
            self._detJ[ke] = det(self.J(*ke))
            return self._detJ[ke]

    def plot(self, **kwds):
        ''' Plots geometry '''
        fig = plt.figure()
        ax = fig.add_axes()
        plt.grid()
        cn = self.corner_nodes
        k = concatenate([linspace(cn[cni][0],cn[(cni+1)%len(cn)][0]) for cni in xrange(len(cn))])
        e = concatenate([linspace(cn[cni][1],cn[(cni+1)%len(cn)][1]) for cni in xrange(len(cn))])
        plt.plot(self.x(k,e),self.y(k,e))
        plt.scatter(self.xcoords,self.ycoords,marker='o',c='b',s=10)
        for (xi,yi,i) in zip(self.xcoords,self.ycoords,range(len(self.xcoords))):
            plt.annotate('%d'%(i+1), xy=(xi,yi), xytext=(0,10), 
            textcoords='offset points', ha='right', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='blue', alpha=0.2))
        if kwds.has_key('filename'):
            plt.savefig(**kwds)
        return fig,ax

class Poi1(Net1D):
    ''' 1-node point. '''
    def __init__(self, *args, **kwds):
        super(Poi1, self).__init__(*args, **kwds)

class Seg2(Net1D):
    ''' 2-node line. '''
    def __init__(self, *args, **kwds):
        super(Seg2, self).__init__(*args, **kwds)
        self.N = lambda k,e: array([0.5*(1-k), 0.5*(1+k)])
        self.dNdk = lambda k,e: array([-0.5,0.5])
        
class Seg3(Net1D):
    ''' 3-node second order line (2 nodes associated with the vertices and 1 with the edge). '''    
    def __init__(self, *args, **kwds):
        super(Seg3, self).__init__(*args, **kwds)
        self.N = lambda ksi,eta: array([ksi*(ksi - 1)/2, ksi*(ksi + 1)/2, -ksi**2 + 1])
        self.dNdk = lambda ksi,eta: array([ksi - 1/2, ksi + 1/2, -2*ksi])

class Tria3(Net2D):
    r""" 3-node triangle. 
        .. image:: img/tria3.png

    SHAPE FUNCTIONS:

    INPUT:
            
    - ``name`` - Name of shape
    - ``[nodes]`` - List of nodes (len(nodes)==3)
        
    OUTPUT:
        
    - Tria3 object
            
    EXAMPLES:

    >>> n1 = Node(x=0.0, y=0.0)
    >>> n2 = Node(x=2.0, y=1.0)
    >>> n3 = Node(x=0.5, y=2.0)
    >>> tr3 = Tria3(name = 'tr3', nodes = [n1, n2, n3])
    
    """
    def __init__(self, *args, **kwds):
        super(Tria3, self).__init__(*args, **kwds)
        self.kecoords = ((0.0,0.0),(1.0,0.0),(0.0,1.0))
        self.corner_nodes = ((0.0,0.0),(1.0,0.0),(0.0,1.0))
        self.N = lambda k,e: array([1-k-e,k,e])
        self.dNdk = lambda k,e: array([-1,1,0])
        self.dNde = lambda k,e: array([-1,0,1])

class Tria6(Net2D):
    r""" 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges. 

    .. image:: img/tria6.png

    SHAPE FUNCTIONS:

    .. math::
        \left\{ \begin{aligned}N_{1}\left(\xi,\eta\right) & =-\left(1-\xi-\eta\right)\left(1-2\left(1-\xi-\eta\right)\right)\\
        N_{2}\left(\xi,\eta\right) & =-\xi\left(1-2\xi\right)\\
        N_{3}\left(\xi,\eta\right) & =-\eta\left(1-2\eta\right)\\
        N_{4}\left(\xi,\eta\right) & =4\xi\left(1-\xi-\eta\right)\\
        N_{5}\left(\xi,\eta\right) & =4\xi\eta\\
        N_{6}\left(\xi,\eta\right) & =4\eta\left(1-\xi-\eta\right)
        \end{aligned}
        \right.

    .. math::
        \left\{ \begin{aligned}\frac{\partial N_{1}\left(\xi,\eta\right)}{\partial\xi} & =1-4\left(1-\xi-\eta\right)\\
        \frac{\partial N_{2}\left(\xi,\eta\right)}{\partial\xi} & =-1+4\xi\\
        \frac{\partial N_{3}\left(\xi,\eta\right)}{\partial\xi} & =0\\
        \frac{\partial N_{4}\left(\xi,\eta\right)}{\partial\xi} & =4\left(1-2\xi-\eta\right)\\
        \frac{\partial N_{5}\left(\xi,\eta\right)}{\partial\xi} & =4\eta\\
        \frac{\partial N_{6}\left(\xi,\eta\right)}{\partial\xi} & =-4\eta
        \end{aligned}
        \right.

    .. math::
        \left\{ \begin{aligned}\frac{\partial N_{1}\left(\xi,\eta\right)}{\partial\eta} & =1-4\left(1-\xi-\eta\right)\\
        \frac{\partial N_{2}\left(\xi,\eta\right)}{\partial\eta} & =0\\
        \frac{\partial N_{3}\left(\xi,\eta\right)}{\partial\eta} & =-1+4\eta\\
        \frac{\partial N_{4}\left(\xi,\eta\right)}{\partial\eta} & =-4\xi\\
        \frac{\partial N_{5}\left(\xi,\eta\right)}{\partial\eta} & =4\xi\\
        \frac{\partial N_{6}\left(\xi,\eta\right)}{\partial\eta} & =4\left(1-\xi-2\eta\right)
        \end{aligned}
        \right.

    INPUT:
            
    - ``name`` - Name of shape
    - ``[nodes]`` - List of nodes (len(nodes)==6)
        
    OUTPUT:
        
    - Tria6 object
            
    EXAMPLES:

    >>> n1 = Node(x=0.0, y=0.0)
    >>> n2 = Node(x=2.0, y=1.0)
    >>> n3 = Node(x=0.5, y=2.0)
    >>> n4 = Node(x=1.0, y=0.0)
    >>> n5 = Node(x=1.0, y=1.0)
    >>> n6 = Node(x=0.2, y=1.0)
    >>> tr6 = Tria6(name = 'tr6', nodes = [n1, n2, n3, n4, n5, n6])
    
    """
    
    def __init__(self, *args, **kwds):
        super(Tria6, self).__init__(*args, **kwds)
        self.kecoords = ((0.0,0.0),(1.0,0.0),(0.0,1.0),(0.5,0.0),(0.5,0.5),(0.0,0.5))
        self.corner_nodes = ((0.0,0.0),(1.0,0.0),(0.0,1.0))
        self.N = lambda ksi, eta: array([
        -(1-ksi-eta)*(1-2*(1-ksi-eta)),
        -ksi*(1-2*ksi),
        -eta*(1-2*eta),
        4*ksi*(1-ksi-eta),
        4*ksi*eta,
        4*eta*(1-ksi-eta),
        ])
        self.dNdk = lambda ksi, eta: array([
        1-4*(1-ksi-eta),
        -1+4*ksi,
        0,
        4*(1-2*ksi-eta),
        4*eta,
        -4*eta,
        ])
        self.dNde = lambda ksi, eta: array([
        1-4*(1-ksi-eta),
        0,
        -1+4*eta,
        -4*ksi,
        4*ksi,
        4*(1-ksi-2*eta),
        ])
    
    def __repr__(self):
        return('{0} : {1}'.format(self.__class__.__name__, self.nodes))


class Quad4(Net2D):
    """ 4-node quadrangle shape
    
    .. image:: img/quad4.png    
    
    SHAPE FUNCTION:

    INPUT:

    OUTPUT:

    EXAMPLES:
    """
    
    def __init__(self, *args, **kwds):
        super(Quad4, self).__init__(*args, **kwds)
        self.kecoords = ((-1,-1),(1,-1),(1,1),(-1,1))
        self.corner_nodes = ((-1,-1),(1,-1),(1,1),(-1,1))
        self.N = lambda ksi, eta: array([
        (1-ksi)*(1-eta)/4,
        (1+ksi)*(1-eta)/4,
        (1+ksi)*(1+eta)/4,
        (1-ksi)*(1+eta)/4,
        ])
        self.dNdk = lambda ksi, eta: array([
        -(1-eta)/4,
        (1-eta)/4,
        (1+eta)/4,
        -(1+eta)/4,
        ])
        self.dNde = lambda ksi, eta: array([
        -(1-ksi)/4,
        -(1+ksi)/4,
        (1+ksi)/4,
        (1-ksi)/4,
        ])

class Quad8(Net2D):
    ''' 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).  '''    
    def __init__(self, *args, **kwds):
        super(Quad8, self).__init__(*args, **kwds)
        self.N = lambda ksi, eta: array([
        (1-ksi)*(1-eta)*(-1-ksi-eta)/4,
        (1+ksi)*(1-eta)*(-1+ksi-eta)/4,
        (1+ksi)*(1+eta)*(-1+ksi+eta)/4,
        (1-ksi)*(1+eta)*(-1-ksi+eta)/4,
        (1-ksi**2)*(1-eta)/2,
        (1+ksi)*(1-eta**2)/2,
        (1-ksi**2)*(1+eta)/2,
        (1-ksi)*(1-eta**2)/2,
        ])
        self.dNdk = lambda ksi, eta: array([
        (1-eta)*(2*ksi+eta)/4,
        (1-eta)*(2*ksi-eta)/4,
        (1+eta)*(2*ksi+eta)/4,
        -(1+eta)*(-2*ksi+eta)/4,
        -ksi*(1-eta),
        (1-eta**2)/2,
        -ksi*(1+eta),
        -(1-eta**2)/2,
        ])
        self.dNde = lambda ksi, eta: array([
        (1-ksi)*(ksi+2*eta)/4,
        -(1+ksi)*(ksi-2*eta)/4,
        (1+ksi)*(ksi+2*eta)/4,
        (1-ksi)*(-ksi+2*eta)/4,
        -(1-ksi**2)/2,
        -eta*(1+ksi),
        (1-ksi**2)/2,
        -eta*(1-ksi),
        ])
        
    def plot(self):
        pass

class Quad9(Net2D):
    ''' 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face). '''    
    def __init__(self, *args, **kwds):
        super(Quad9, self).__init__(*args, **kwds)
        self.N = lambda ksi, eta: array([
        ksi*eta*(ksi-1)*(eta-1)/4,
        ksi*eta*(ksi+1)*(eta-1)/4,
        ksi*eta*(ksi+1)*(eta+1)/4,
        ksi*eta*(ksi-1)*(eta+1)/4,
        (1-ksi**2)*eta*(eta-1)/2,
        ksi*(ksi+1)*(1-eta**2)/2,
        (1-ksi**2)*eta*(eta+1)/2,
        ksi*(ksi-1)*(1-eta**2)/2,
        (1-ksi**2)*(1-eta**2),
        ])
        self.dNdk = lambda ksi, eta: array([
        (2*ksi-1)*eta*(eta-1)/4,
        (2*ksi+1)*eta*(eta-1)/4,
        (2*ksi+1)*eta*(eta+1)/4,
        (2*ksi-1)*eta*(eta+1)/4,
        -ksi*eta*(eta-1),
        (2*ksi+1)*(1-eta**2)/2,
        -ksi*eta*(eta+1),
        (2*ksi-1)*(1-eta**2)/2,
        -2*ksi*(1-eta**2),
        ])
        self.dNde = lambda ksi, eta: array([
        ksi*(ksi-1)*(2*eta-1)/4,
        ksi*(ksi+1)*(2*eta-1)/4,
        ksi*(ksi+1)*(2*eta+1)/4,
        ksi*(ksi-1)*(2*eta+1)/4,
        (1-ksi**2)*(2*eta-1)/2,
        -ksi*eta*(ksi+1),
        (1-ksi**2)*(2*eta+1)/2,
        -ksi*eta*(ksi-1),
        -2*eta*(1-ksi**2),
        ])

