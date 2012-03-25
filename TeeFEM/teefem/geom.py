# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 23:12:23 2012

@author: Jukka Aho

Mesh routines

"""

from __future__ import division


import teefem
from teefem import log
import os
import re
from numpy import array, matrix, concatenate, linspace
from numpy.linalg import inv, det

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
        self.groups = set()
        self.boundary_conditions = set()
        self.bc = {}
        self.load = {}
        
    def __repr__(self):
        return('{0} : ({1.x:0.2f},{1.y:0.2f},{1.z:0.2f})'.format(self.__class__.__name__, self))
    
    def assign_boundary_condition(self, bc):
        self.boundary_conditions.add(bc)

    def update_boundary_conditions(self):
        for bc in self.boundary_conditions:
            dim = ('dx','dy','dz','drx','dry','drz')
            for i in range(len(dim)):
                if bc.exist(dim[i]):
                    self.bc[dim[i]] = bc.val(dim[i])
            dim = ('fx','fy','fz','mx','my','mz')
            for i in range(len(dim)):
                if bc.exist(dim[i]):
                    self.load[dim[i]] = bc.val(dim[i])
    
    def hasbc(self, dof):
        return self.bc.has_key(dof)
    
    def valbc(self, dof):
        return self.bc[dof]

    def hasload(self, dof):
        return self.load.has_key(dof)
    
    def valload(self, dof):
        return self.load[dof]
    
    def update_field(self, *args, **kwds):
        field = kwds['field'].upper()
        if not self.fields.has_key(field):
            self.fields[field] = {}
        params = kwds.get('params')
        values = kwds.get('values')
#        if not isinstance(params, list):
#            params = list([params])
#        if not isinstance(values, list):
#            values = list([values])
        for (p,v) in zip(params, values):
            self.fields[field][p] = v
    
    def status(self):
        print self.bc
        print self.load
        print self.boundary_conditions


class Net(object):
    
    count = 0    
    
    def __init__(self, *args, **kwds):
        ''' General Net object
        
        EXAMPLES:        
        
        >>> node1 = Node(x=1,y=2,z=3)
        >>> node2 = Node(x=2,y=3,z=4)
        >>> node3 = Node(x=3,y=4,z=5)
        >>> net1 = Net(name = 'net1', nodes = [node1, node2, node3])
        
        '''
        Net.count += 1
        self.groups = set()
        self.name = kwds.get('name', 'shape%d'%(Net.count))
        self.nodes = kwds.get('nodes')
        self.xcoords = [n.x for n in self.nodes]
        self.ycoords = [n.y for n in self.nodes]
        self.zcoords = [n.z for n in self.nodes]
        self.x = lambda *ke: sum(N*nodes.x for (N,nodes) in zip(self.N(*ke),self.nodes))
        self.y = lambda *ke: sum(N*nodes.y for (N,nodes) in zip(self.N(*ke),self.nodes))
        self.z = lambda *ke: sum(N*nodes.z for (N,nodes) in zip(self.N(*ke),self.nodes))

class Net1D(Net):
    """
    1d shapes
    """
    def __init__(self, *args, **kwds):
        super(Net1D, self).__init__(*args, **kwds)
        self.Jk = lambda k: sum(dNdk*nodes.x for (dNdk,nodes) in zip(self.dNdk(k),self.nodes))
        self.Je = lambda k: sum(dNdk*nodes.y for (dNdk,nodes) in zip(self.dNdk(k),self.nodes))
        self.invJk = lambda k: 1/self.Jk(k)
        self.invJe = lambda k: 1/self.Je(k)
#
#        self.J11 = lambda *ke: sum(dNdk*nodes.x for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))
#        self.J12 = lambda *ke: sum(dNdk*nodes.y for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))

class Net2D(Net):
    
    """ 2d shapes """
    
    def __init__(self, *args, **kwds):
        super(Net2D, self).__init__(*args, **kwds)

    def J(self,*ke):
        """ 
        2d Jacobian matrix 
        """
        J11 = sum(dNdk*nodes.x for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))
        J12 = sum(dNdk*nodes.y for (dNdk,nodes) in zip(self.dNdk(*ke),self.nodes))
        J21 = sum(dNde*nodes.x for (dNde,nodes) in zip(self.dNde(*ke),self.nodes))
        J22 = sum(dNde*nodes.y for (dNde,nodes) in zip(self.dNde(*ke),self.nodes))
        return matrix([[J11,J12],[J21,J22]])
    
    def invJ(self,*ke):
        ''' Returns inverse Jacobian matrix '''
        return inv(self.J(*ke))
    
    def detJ(self,*ke):
        ''' Jacobian determinant '''
        return det(self.J(*ke))

    def plot(self, **kwds):
        ''' Plot geometry '''
        import matplotlib.pylab as plt
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

    def plot3d(self, **kwds):
        ''' Plot geometry with mplot3d '''
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.grid()
        cn = self.corner_nodes
        k = concatenate([linspace(cn[cni][0],cn[(cni+1)%len(cn)][0]) for cni in xrange(len(cn))])
        e = concatenate([linspace(cn[cni][1],cn[(cni+1)%len(cn)][1]) for cni in xrange(len(cn))])
        ax.plot(self.x(k,e),self.y(k,e),self.z(k,e))
        ax.scatter(self.xcoords,self.ycoords,self.zcoords,marker='o',c='b',s=10)
        for (xi,yi,zi,i) in zip(self.xcoords,self.ycoords,self.zcoords,range(len(self.xcoords))):
            label = '%d'%(i+1)
            ax.text(xi,yi,zi,label)
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
        self.N = lambda k: array([0.5*(1-k), 0.5*(1+k)])
        self.dNdk = lambda k: array([-0.5,0.5])
        
class Seg3(Net1D):
    ''' 3-node second order line (2 nodes associated with the vertices and 1 with the edge). '''    
    def __init__(self, *args, **kwds):
        super(Seg3, self).__init__(*args, **kwds)
        self.N = lambda ksi: array([ksi*(ksi - 1)/2, ksi*(ksi + 1)/2, -ksi**2 + 1])
        self.dNdk = lambda ksi: array([ksi - 1/2, ksi + 1/2, -2*ksi])
        

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
        x1 = self.xcoords[0]
        x2 = self.xcoords[1]
        x3 = self.xcoords[2]
        y1 = self.ycoords[0]
        y2 = self.ycoords[1]
        y3 = self.ycoords[2]
        x21 = x2-x1
        y21 = y2-y1
        x31 = x3-x1
        y31 = y3-y1
        self.area = 0.5*(x21*y31 - y21*x31)
        
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


def parsemail(fh):
    
    ''' Parse Aster mail format '''
    
    mesh = Mesh()

    def titre(*args):
        pass

    def coor_2d(*args):
        nbr = 3
        res = int(len(args)/nbr)
        for i in xrange(res):
            mesh.nodes[args[nbr*i]] = Node(x = float(args[nbr*i+1]), y = float(args[nbr*i+2]))

    def coor_3d(*args):
        nbr = 4
        res = int(len(args)/nbr)
        for i in xrange(res):
            mesh.nodes[args[nbr*i]] = Node(x = float(args[nbr*i+1]), y = float(args[nbr*i+2]), z = float(args[nbr*i+3]))

    def init_net(**kwds):
        params = kwds['params']
        items = kwds['items']
        nbr = params['nbr']
        res = int(len(items)/nbr)
        for i in xrange(res):
            netname = items[nbr*i]
            nodenames = items[nbr*i+1:nbr*(i+1)]
            mesh.nets[netname] = params['netf'](name = netname, nodes = [mesh.nodes[ni] for ni in nodenames])

    def group_no(*args):
        name = args[0]
        mesh.group_no[name] = set([mesh.nodes[ni] for ni in args[1:]])
        for ni in args[1:]:
            mesh.nodes[ni].groups.add(name)

    def group_ma(*args):
        name = args[0]
        mesh.group_ma[name] = set([mesh.nets[ni] for ni in args[1:]])
        for ni in args[1:]:
            mesh.nets[ni].groups.add(name)

    def _F(**kwds):
        return kwds

    functions = {
        'TITRE' : titre,
        'COOR_2D': coor_2d,
        'COOR_3D': coor_3d,
        'GROUP_NO': group_no,
        'GROUP_MA': group_ma,
    }

    functions2 = {
        'POI1' : _F(netf=Poi1, nbr=2),
        'SEG2' : _F(netf=Seg2, nbr=3),
        'SEG3' : _F(netf=Seg3, nbr=4),
        'TRIA3' : _F(netf=Tria3, nbr=4),
        'TRIA6' : _F(netf=Tria6, nbr=7),
        'QUAD4' : _F(netf=Quad4, nbr=5),
        'QUAD8' : _F(netf=Quad8, nbr=9),
        'QUAD9' : _F(netf=Quad9, nbr=10),
    }
    
    data = fh.read()
    for section in re.split('FINSF', data):
        items = re.split('\s+', section)
        items = [i.strip('% ') for i in items]
        items = [i for i in items if len(i) != 0]
        if items[0] in functions:
            functions[items[0]](*items[1:])
        if items[0] in functions2:
            init_net(items = items[1:], netname = items[0], params = functions2[items[0]])
    
    return mesh

def parsemsh(fh):
    
    ''' Read Gmsh .msh file format '''
    
    mesh = Mesh()

    mapping = {
        1: Seg2, 
        2: Tria3, 
        3: Quad4,
        8: Seg3,
        9: Tria6,
        10: Quad9,
        15: Poi1,
        16: Quad8,
        }
    regions = {}
    section = None
    for line in fh:
        line = line.strip()
        if line.startswith('$'):
            if line.startswith('$End'): section = None
            else: section = line
            continue
        columns = line.split()
        if len(columns)==1: continue # Emme tarvitse lukumäärää                
        if section == "$MeshFormat": mesh.meshformat = line
        if section == "$PhysicalNames":
            name = columns[2].strip('"')
            mesh.group_ma[name] = set()
#            regions[int(columns[1])] = mesh.group_ma[name]
            regions[int(columns[1])] = name
        if section == "$Nodes":
            mesh.nodes['N%d'%(int(columns[0]))] = Node(x = float(columns[1]),y = float(columns[2]),z = float(columns[3]))
        if section == "$Elements":
            columns = [int(ci) for ci in columns]
            eid,etypeid,enumofparams,eregion,egentity = columns[0:5]
#                    eparams = columns[4:enumofparams-2]
            enodes = columns[enumofparams+3:]
#                    self.elements[eid] = mapping[etypeid](self.regions[eregion],egentity,eparams,[self.nodes[ni] for ni in enodes])
            try:
                shape = mapping[etypeid](name = 'M%d'%(int(eid)), nodes = [mesh.nodes['N%d'%(ni)] for ni in enodes])
            except KeyError:
                print("Unable to create shape type {0}".format(etypeid))
                print columns
            mesh.nets['M%d'%(int(eid))] = shape
            shape.groups.add(regions[eregion])
            mesh.group_ma[regions[eregion]].add(shape)
            
    return mesh
#    mesh.create_node_groups()


def parsemed(fh):
    ''' Read Salome .med file format '''
    print("I'm broken")
    return Mesh()
    
    db = tables.openFile(kwds['filename'])
    NOE = db.getNode('/ENS_MAA/geom/NOE')
    MAI = db.getNode('/ENS_MAA/geom/MAI')
    FAS = db.getNode('/ENS_MAA/geom/FAS')
    # GROUP_MA
    for group_ma in FAS.ELEME._f_iterNodes():
        name = ''.join(group_ma.GRO.NOM[:]).strip()
        idx = group_ma._v_attrs.NUM
        GMA = GROUP_MA(name, idx)
        console_print("New group_ma: {0} with index {1}".format(name,idx))
        self.group_ma_idx[idx] = GMA
        self.group_ma_name[name] = GMA
    # NODES
    COO = NOE.COO[:]
    FAM = NOE.FAM[:]
    NBR = NOE.COO._v_attrs.NBR
    print("Number of nodes: %d"%(NBR))
    for idx in xrange(NBR):
        fam = FAM[idx]
        coo = COO[idx::NBR]
        node = Node(coo)
        self.nodes[idx+1] = node
    # NETS
    mapping = {
        "PO1": POI1, 
        "SE2": SEG2,
        "SE3": SEG3,
        "TR3": TRIA3,
        "TR6": TRIA6,
        }
    for (typ,arrs) in MAI._v_children.iteritems():
        NBR = int(arrs.NOD._v_attrs.NBR)
        console_print("Number of {type} nets: {nbr:d}".format(type=typ, nbr=NBR))
        ntyp = mapping[typ]
        NOD = arrs.NOD[:]
        FAM = arrs.FAM[:]
        for idx in xrange(NBR):
            fam = FAM[idx]
            nod = NOD[idx::NBR]
#                net = ntyp(nodes = [self.nodes[ni] for ni in nod][::-1], group_ma = self.group_ma_idx[fam])
            net = ntyp(nodes = [self.nodes[ni] for ni in nod], group_ma = self.group_ma_idx[fam])
            self.nets.add(net)


class Mesh(object):
    
    ''' Mesh class. Reads file formats: Gmsh .msh, Salome .med, Aster .mail '''    
    
    def __init__(self, *args, **kwds):
        ''' Initializes mesh class '''
        self.group_no = {}
        self.group_ma = {}
        self.nodes = {}
        self.nets = {}

    def find_group_ma(self, gma):
        return self.group_ma[gma]

    def find_group_no(self, gno):
        return self.group_ma[gno]
        
    def create_node_groups(self):
        ''' Create node groups from element groups '''
        for key in self.group_ma.iterkeys():
            self.group_no[key] = set()
            for element in self.group_ma[key]:
                for node in element.nodes:
                    node.groups.add(key)
                    self.group_no[key].add(node)
    
    @teefem.deprecated
    def assign_material(self, mat):
        for n in self.nets.itervalues():
            n.material = mat
    
    def status(self, fulloutput=False):
        ret = ''
        ret += 'Shape groups:\n'
        for group_ma in self.group_ma:
            ret += '{0:10s} : {1} shapes\n'.format(group_ma, len(self.group_ma[group_ma]))
        ret += 'Node groups:\n'
        for group_no in self.group_no:
            ret += '{0:10s} : {1} nodes\n'.format(group_no, len(self.group_no[group_no]))
        if fulloutput:
            ret += 'Nodes: \n'
            for (key, val) in self.nodes.iteritems():
                ret += '{0} : (x,y,z) -> ({1},{2},{3})\n'.format(key, val.x, val.y, val.z)
        return ret

def mesh(**kwds):
    fileformat = os.path.splitext(kwds['filename'])[-1]
    func = {
        '.msh': parsemsh,
#        '.med': parsemed,
        '.mail': parsemail,
        }
    with open(kwds['filename'], 'r') as filehandler:
        return func[fileformat](filehandler)

log.info("Module {0} loaded.".format(__file__))
