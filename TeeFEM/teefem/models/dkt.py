# -*- coding: utf-8 -*-
"""
===========================
Discrete Kirchhoff elements
===========================

Examples
--------
- `examples/dkt_ex1.py` -- dkt_ex1.py

"""

from __future__ import division

import teefem
import teefem.elements
from teefem.boundary_conditions import PressureBoundaryCondition
from teefem.plate_functions import PlateCharacteristic
from common import Model
from teefem import log, cache

#from scipy.sparse.linalg import spsolve
from numpy import array, zeros, concatenate, linspace, matrix
from numpy import sqrt
import numpy as np
#import matplotlib.tri as tri
np.set_printoptions(precision = 2, edgeitems=20)

def checkmatrix(A):
    teefem.log.info('Determinant: {0}'.format(np.linalg.det(A)))
    s = A.shape
    teefem.log.info('Shape: {0}'.format(s))
    for i in range(s[0]):
        print("Sum col {0}: {1}".format(i+1, np.sum(A[i,:])))

class MEBODKT(teefem.elements.Element1D):
    ''' SEG2 element in DKT formulation '''
    def __init__(self, *args, **kwds):
        self.dimension = 2*3
        super(MEBODKT, self).__init__(*args, **kwds)

###############################################################################

class MEDKTR3(teefem.elements.Element2D):
    
    ''' DKT finite element '''
    
    def __init__(self, *args, **kwds):
        super(MEDKTR3, self).__init__(*args, **kwds)
        
        self.boundary_conditions = set()
        
        self.dimension = 3*3
        self.degrees_of_freedom = ('DZ','DRX','DRY')
        self.geom = kwds.get('geom')
        # self.ipoints = [(1.0/3.0, 1.0/3.0)]
        # self.iweights = [(1.0/2.0)]
        self.ipoints = [(0.5, 0.0), (0.0, 0.5), (0.5, 0.5)]
        self.iweights = [(1.0/6.0), (1.0/6.0), (1.0/6.0)]
        
        # Shape functions
        
        self.Lxy = lambda k,e: array([1-k-e,k,e])
        
        N1 = lambda k,e: 2*(1-k-e)*(0.5*k-e)
        N2 = lambda k,e: k*(2*k-1)
        N3 = lambda k,e: e*(2*e-1)
        N4 = lambda k,e: 4*k*e
        N5 = lambda k,e: 4*e*(1-k-e)
        N6 = lambda k,e: 4*k*(1-k-e)
        
        x1 = self.geom.nodes[0].x
        x2 = self.geom.nodes[1].x
        x3 = self.geom.nodes[2].x
        y1 = self.geom.nodes[0].y
        y2 = self.geom.nodes[1].y
        y3 = self.geom.nodes[2].y
        x12 = x1-x2
        x31 = x3-x1
        x23 = x2-x3
        y12 = y1-y2
        y31 = y3-y1
        y23 = y2-y3
        l23 = x23**2 + y23**2
        l31 = x31**2 + y31**2
        l12 = x12**2 + y12**2

        a4 = -x23/l23
        b4 = 3/4*x23*y23/l23
        c4 = (1/4*x23**2-1/2*y23**2)/l23
        d4 = -y23/l23
        e4 = (1/4*y23**2 - 1/2*x23**2)/l23

        a5 = -x31/l31
        b5 = 3/4*x31*y31/l31
        c5 = (1/4*x31**2-1/2*y31**2)/l31
        d5 = -y31/l31
        e5 = (1/4*y31**2 - 1/2*x31**2)/l31

        a6 = -x12/l12
        b6 = 3/4*x12*y12/l12
        c6 = (1/4*x12**2-1/2*y12**2)/l12
        d6 = -y12/l12
        e6 = (1/4*y12**2 - 1/2*x12**2)/l12
        
        Hx1 = lambda k,e: 1.5*(a6*N6(k,e)-a5*N5(k,e))
        Hx2 = lambda k,e: b5*N5(k,e) + b6*N6(k,e)
        Hx3 = lambda k,e: N1(k,e) - c5*N5(k,e) - c6*N6(k,e)
        Hy1 = lambda k,e: 1.5*(d6*N6(k,e) - d5*N5(k,e))
        Hy2 = lambda k,e: -N1(k,e) + e5*N5(k,e) + e6*N6(k,e)
        Hy3 = lambda k,e: -b5*N5(k,e) - b6*N6(k,e)

        Hx4 = lambda k,e: 1.5*(a4*N4(k,e)-a6*N6(k,e))
        Hx5 = lambda k,e: b6*N6(k,e) + b4*N4(k,e)
        Hx6 = lambda k,e: N2(k,e) - c6*N6(k,e) - c4*N4(k,e)
        Hy4 = lambda k,e: 1.5*(d4*N4(k,e) - d6*N6(k,e))
        Hy5 = lambda k,e: -N2(k,e) + e6*N6(k,e) + e4*N4(k,e)
        Hy6 = lambda k,e: -b6*N6(k,e) - b4*N4(k,e)

        Hx7 = lambda k,e: 1.5*(a5*N5(k,e)-a4*N4(k,e))
        Hx8 = lambda k,e: b4*N4(k,e) + b5*N5(k,e)
        Hx9 = lambda k,e: N3(k,e) - c4*N4(k,e) - c5*N5(k,e)
        Hy7 = lambda k,e: 1.5*(d5*N5(k,e) - d4*N4(k,e))
        Hy8 = lambda k,e: -N3(k,e) + e4*N4(k,e) + e5*N5(k,e)
        Hy9 = lambda k,e: -b4*N4(k,e) - b5*N5(k,e)
        
        self.Hx = lambda k,e: array([Hx1(k,e), Hx2(k,e), Hx3(k,e), Hx4(k,e), Hx5(k,e), Hx6(k,e), Hx7(k,e), Hx8(k,e), Hx9(k,e)])
        self.Hy = lambda k,e: array([Hy1(k,e), Hy2(k,e), Hy3(k,e), Hy4(k,e), Hy5(k,e), Hy6(k,e), Hy7(k,e), Hy8(k,e), Hy9(k,e)])
        
        self.Bx = lambda k,e: np.sum(self.Hx(k,e)*self.U)
        self.By = lambda k,e: np.sum(self.Hy(k,e)*self.U)

    def assign_material(self, mat):
        self.material = mat
#        log.debug("Material assigned. E: {0}".format(self.material.E))

    def assign_boundary_condition(self, bc):
        if bc.__class__ is PressureBoundaryCondition:
            self.pressure = bc.pressure

    def assign_char(self, char):
        if char.__class__ is PlateCharacteristic:
            self.thickness = char.thickness

    @property
    def U(self):
        return array([[
            n.fields['DEPL']['DZ'], 
            n.fields['DEPL']['DRX'], 
            n.fields['DEPL']['DRY']] for n in self.geom.nodes]).flatten()
    
    def w(self,k,e):
        ''' Jostakin löytyi tämmöset muofofunktioit joilla voisi yrittää
        interpoloida siirtymäkenttää elementin alueella. '''
        U = self.U
        la = 1-k-e
        N1 = 3*la**2 - 2*la**3 + 2*k*e*la
        N2 = la**2*k + k*e*la/2
        N3 = la**2*e + k*e*la/2
        N4 = 3*k**2 - 2*k**3 + 2*k*e*la
        N5 = k**2*(k-1) - k*e*la
        N6 = k**2*e + k*e*la/2
        N7 = 3*e**2 - 2*e**3 + 2*k*e*la
        N8 = e**2*k + k*e*la/2
        N9 = e**2*(e-1) - k*e*la
        N = array([N1,N2,N3,N4,N5,N6,N7,N8,N9])
        return np.sum(N*U)

    @property
    @cache
    def detJ(self):
        return self.geom.detJ

    @cache
    def B(self,*ke):
        k = ke[0]
        e = ke[1]

        x1 = self.geom.nodes[0].x
        x2 = self.geom.nodes[1].x
        x3 = self.geom.nodes[2].x
        y1 = self.geom.nodes[0].y
        y2 = self.geom.nodes[1].y
        y3 = self.geom.nodes[2].y

        x12 = x1-x2
        x31 = x3-x1
        x23 = x2-x3
        y12 = y1-y2
        y31 = y3-y1
        y23 = y2-y3

        A2 = x31*y12-x12*y31
        self.A = A2/2
        l12 = sqrt(x12**2+y12**2)
        l23 = sqrt(x23**2+y23**2)
        l31 = sqrt(x31**2+y31**2)
        P4 = -6*x23/l23**2
        P5 = -6*x31/l31**2
        P6 = -6*x12/l12**2
        t4 = -6*y23/l23**2
        t5 = -6*y31/l31**2
        t6 = -6*y12/l12**2
        q4 = 3*x23*y23/l23**2
        q5 = 3*x31*y31/l31**2
        q6 = 3*x12*y12/l12**2
        r4 = 3*y23**2/l23**2
        r5 = 3*y31**2/l31**2
        r6 = 3*y12**2/l12**2
        
        Hxk = array([
            P6*(1-2*k)+(P5-P6)*e, 
            q6*(1-2*k)-(q5+q6)*e, 
            -4+6*(k+e)+r6*(1-2*k)-e*(r5+r6),
            -P6*(1-2*k)+e*(P4+P6),
            q6*(1-2*k)-e*(q6-q4),
            -2+6*k+r6*(1-2*k)+e*(r4-r6),
            -e*(P5+P4),
            e*(q4-q5),
            -e*(r5-r4)])
            
        Hyk = array([
            t6*(1-2*k)+e*(t5-t6),
            1+r6*(1-2*k)-e*(r5+r6),
            -q6*(1-2*k)+e*(q5+q6),
            -t6*(1-2*k)+e*(t4+t6),
            -1+r6*(1-2*k)+e*(r4-r6),
            -q6*(1-2*k)-e*(q4-q6),
            -e*(t4+t5),
             e*(r4-r5),
            -e*(q4-q5)])

        Hxe = array([
            -P5*(1-2*e)-k*(P6-P5),
              q5*(1-2*e)-k*(q5+q6),
             -4+6*(k+e)+r5*(1-2*e)-k*(r5+r6),
              k*(P4+P6),
              k*(q4-q6),
             -k*(r6-r4),
              P5*(1-2*e)-k*(P4+P5),
              q5*(1-2*e)+k*(q4-q5),
             -2+6*e+r5*(1-2*e)+k*(r4-r5)])

        Hye = array([
            -t5*(1-2*e)-k*(t6-t5),
            1+r5*(1-2*e)-k*(r5+r6),
            -q5*(1-2*e)+k*(q5+q6),
            k*(t4+t6),
            k*(r4-r6),
            -k*(q4-q6),
            t5*(1-2*e)-k*(t4+t5),
            -1+r5*(1-2*e)+k*(r4-r5),
            -q5*(1-2*e)-k*(q4-q5)])

        return (1/A2)*array([
            y31*Hxk+y12*Hxe, 
            -x31*Hyk-x12*Hye, 
            -x31*Hxk-x12*Hxe+y31*Hyk+y12*Hye])

    @cache
    def D(self, *ke):
        try:
            material = self.material
        except AttributeError as err:
            log.error(err)
            log.error('Element groups: {0}'.format(self.groups))
        E = material.E
        nu = material.nu
        t = self.thickness(*ke)
        return E*t**3/(12*(1-nu**2)) * matrix([[1,nu,0],[nu,1,0],[0,0,0.5*(1-nu)]])

    @property
    def stiffness_matrix(self):
        ''' Stiffness matrix '''
        K = matrix(zeros((self.dimension,self.dimension)))
        for (W,ke) in zip(self.iweights, self.ipoints):
            detJ = self.detJ(*ke)
            B = self.B(*ke)
            D = self.D(*ke)
            K += W*B.T*D*B*detJ
        return K

    @property
    def force_vector(self):
        ''' Palauttaa kuormitusvektorin '''
        b = array([W*self.pressure(*ke)*self.Lxy(*ke)*self.detJ(*ke) for (W,ke) in zip(self.iweights, self.ipoints)])
        F = zeros(self.dimension)
        F[0::3] = np.sum(b,axis=0)
        return matrix(F).T

    def plot(self, **kwds):
        ''' Plot geometry with mplot3d '''
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.grid()
        cn = self.geom.corner_nodes
        k = concatenate([linspace(cn[cni][0],cn[(cni+1)%len(cn)][0]) for cni in xrange(len(cn))])
        e = concatenate([linspace(cn[cni][1],cn[(cni+1)%len(cn)][1]) for cni in xrange(len(cn))])
        w = np.vectorize(self.w)
        ax.plot(self.geom.x(k,e),self.geom.y(k,e),self.geom.z(k,e)+w(k,e))
        ax.scatter(self.geom.xcoords,self.geom.ycoords,self.geom.zcoords,marker='o',c='b',s=10)
        for (xi,yi,zi,i) in zip(self.geom.xcoords,self.geom.ycoords,self.geom.zcoords,range(len(self.geom.xcoords))):
            label = '%d'%(i+1)
            ax.text(xi,yi,zi,label)
        if kwds.has_key('filename'):
            plt.savefig(**kwds)
        return fig,ax    

    def plot2(self, **kwds):
        ''' Plot geometry with mplot3d '''
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.grid()

        X = np.linspace(0,1,25)
        Xran = range(len(X))
        w = lambda k,e: 0*k*np.sin(e)
        data = np.array([(X[i],X[j],self.w(X[i],X[j])) for j in Xran for i in Xran[0:len(X)-j]])
        tr = tri.Triangulation(data[:,0], data[:,1])
        for edge in tr.edges:
            idx1 = edge[0]
            idx2 = edge[1]
            x1 = data[idx1,0]
            y1 = data[idx1,1]
            z1 = data[idx1,2]
            x2 = data[idx2,0]
            y2 = data[idx2,1]
            z2 = data[idx2,2]
            ax.plot([x1,x2],[y1,y2],[z1,z2], 'b')
        return fig,ax    
        
    def update(self, U):
        ''' '''
        u = zeros(self.dimension)
#        print("dim: {0}".format(u.shape))
        for i in xrange(len(self.nodes)):
            node = self.nodes[i]
            node.update_field(
                field = 'DEPL', 
                params = self.degrees_of_freedom, 
                values = (U[node.gdof[0]], U[node.gdof[1]], U[node.gdof[2]]))
            for j in xrange(3):
                u[i*3+j] = U[node.gdof[j]]
        self.u = matrix(u).T

    @cache
    def M(self,*ke):
        ''' Momentti pisteessä (k,e) '''
        B = self.B(*ke)
        D = self.D(*ke)
        return D*B*self.u

    @cache
    def S(self, *ke):
        ''' Jännitys pisteessä (k,e) '''
        return 6/self.thickness(*ke)**2*self.M(*ke)

    @cache
    def vmis(self, *ke):
        S = self.S(*ke)
#        print S
        sxx = S[0,0]
        syy = S[1,0]
        sxy = S[2,0]
        return np.sqrt(sxx**2 + syy**2 - sxx*syy + 3*sxy**2)

###############################################################################

class DKT(Model):
    
    ''' Discrete Kirchhoff Triangle model. '''
    
    def __init__(self, *args, **kwds):
        super(DKT, self).__init__(*args, **kwds)
        self.mapping = {
            'Seg2': MEBODKT,
            'Tria3': MEDKTR3,
        }
        self.nodedim = 3
        self.nodedofs = ('dz','drx','dry')
        self.nodeloads = ('fz','mx','my')
        self.init()

    def plot(self, **kwds):
        ''' Plottaus. '''
        import matplotlib.pylab as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        ax = fig.gca(projection = '3d')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        ax.grid()
#        x = []
#        y = []
#        z = []
#        for n in self.nodes:
#            x.append(n.x)
#            y.append(n.y)
#            z.append(n.z+n.fields['DEPL']['DZ'])
#        ax.scatter(x,y,z,marker='o',c='b',s=10)
        for e in self.stiff_elements:
            x1 = e.nodes[0].x
            y1 = e.nodes[0].y
            z1 = e.nodes[0].z + e.nodes[0].fields['DEPL']['DZ']
            x2 = e.nodes[1].x
            y2 = e.nodes[1].y
            z2 = e.nodes[1].z + e.nodes[1].fields['DEPL']['DZ']
            x3 = e.nodes[2].x
            y3 = e.nodes[2].y
            z3 = e.nodes[2].z + e.nodes[2].fields['DEPL']['DZ']
            ax.plot([x1,x2,x3,x1],[y1,y2,y3,y1],[z1,z2,z3,z1],'--ko')
        return fig,ax
        
dkt = DKT

###############################################################################
###############################################################################
###############################################################################

#def ex2():
#    
#    ''' DKT esimerkki 2, MEDKTR3-elementin formuloinnin tarkastus '''
#    
#    n1 = teefem.geom.Node(x=0, y=0, z=0.0)
#    n2 = teefem.geom.Node(x=1, y=0, z=0.0)
#    n3 = teefem.geom.Node(x=0.5, y=0.5, z=0.0)
#    g1 = teefem.geom.Tria3(nodes = (n1,n2,n3))
#    e1 = MEDKTR3(geom = g1)
#    print e1.status
#    import matplotlib.pylab as plt
##    g1.plot3d()
#    # Aiheutaan jokin satunnainen siirtymätila
#    n1.update_field(field = 'DEPL', params = ('DZ','DRX','DRY'), values = (0.0,0.0,1.0))
#    n2.update_field(field = 'DEPL', params = ('DZ','DRX','DRY'), values = (0.0,1.0,0.0))
#    n3.update_field(field = 'DEPL', params = ('DZ','DRX','DRY'), values = (1.0,0.0,0.0))
#    print n1.fields
#    print e1.U
#    print e1.Hx(1/2,1/2)
#    print e1.Bx(1/2,1/2)
#    print e1.By(1/2,1/2)
#    print("Disp: {0}".format(e1.w(1,0)))
#    #e1.plot2()
#    #plt.show()
#    print e1.kinematic_matrix(1/2,1/2)
#    g1.material = teefem.materials.Elastic(E = 100.0e9, nu = 0.5)
#    print e1.material_matrix
#    k1 = e1.stiffness_matrix/1e8
#    print k1
#    # print k1.T - k1
#    # print np.linalg.det(k1)

def plotmdl(mdl, vmis):
    ''' Plottailee kaikenlaisia kuvaajia '''

    # Mises jännitys integroimispisteitä interpoloimalla
    # http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata

#    plt.figure(1)
#    mdl.plot()
    
    fig2 = plt.figure()
    x = []
    y = []
    Mx = []
    My = []
    Mxy = []
    Sx = []
    Sy = []
    Sxy = []
    Vmis = []
    Vmis2 = []
    for e in mdl.elset['OM1']:
        for ke in e.ipoints:
            x.append(e.geom.x(*ke))
            y.append(e.geom.y(*ke))
            M = e.M(*ke)
            Mx.append(M[0,0])
            My.append(M[1,0])
            Mxy.append(M[2,0])
            S = e.S(*ke)/1e6
            Sx.append(S[0,0])
            Sy.append(S[1,0])
            Sxy.append(S[1,0])
            Vmis.append(e.vmis(*ke)/1e6)

    print("Mises max: {0}".format(max(Vmis)))
    
    for i in range(len(x)):
        r = np.sqrt(x[i]**2 + y[i]**2)
        vm = vmis(r)/1e6
        Vmis2.append(vm)
        print("r: {0:0.2f} Vmis(FEM): {1:0.2f} Vmis(Acc): {2:0.2f}".format(r,Vmis[i],vm))

    Vmis_FEM = np.array(Vmis)
    Vmis_ACC = np.array(Vmis2)
    

    # define grid.
    xi = np.linspace(min(x),max(x),100)
    yi = np.linspace(min(y),max(y),100)
    # grid the data.
    plt.subplot(121)
    si = griddata((x, y), Vmis, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,si,30,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,si,30,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
    # plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(min(x),max(x))
    plt.ylim(min(y),max(y))
    plt.title('FEM')
    #plt.tight_layout()

    plt.subplot(122)
    si = griddata((x, y), Vmis2, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,si,30,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,si,30,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
    # plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(min(x),max(x))
    plt.ylim(min(y),max(y))
    plt.title('Tarkka')

    plt.show()


def ex1():
    ''' Suorakaidelaatta, tasan jakautunut kuorma '''
    mesh = teefem.mesh.unitsquare()
    print mesh.status(fulloutput = True)
    mdl = DKT(mesh = mesh)
    mat = teefem.materials.elastic(E = 100.0e9, nu = 0.3)
    teefem.assign_material(elements = mdl.elset['OM1'], material = mat)
    load = teefem.pressure_bc(pressure = lambda k,e: -100.0e3)
    encastre = teefem.dirichlet_bc(encastre = True)
    teefem.assign_bc(elements = mdl.elset['OM1'], bc = load)
    teefem.assign_bc(nodes = mdl.nset['GA1'], bc = encastre)
    teefem.assign_bc(nodes = mdl.nset['GA2'], bc = encastre)
    teefem.assign_bc(nodes = mdl.nset['GA3'], bc = encastre)
    teefem.assign_bc(nodes = mdl.nset['GA4'], bc = encastre)
    carel = teefem.plate_functions.PlateCharacteristic(thickness = lambda k,e: 20e-3)
    teefem.assign_char(elements = mdl.elset['OM1'], char = carel)
    mdl.static_solve()
    dymin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])
    print("DZ: %0.14E"%(dymin*1000))
    plotmdl(mdl)

def rotsy():
    ''' Ympyrälaatan tarkkoja ratkaisuja '''
    import sympy
    q,a,r,E,t,nu,D = sympy.var('q,a,r,E,t,nu,D')
#    D = E*t**3/(12*(1-nu**2))
    w = -q*a**4/(64*D)*(1-r**2/a**2)*(1-r**2/a**2)
    print("w: {0}".format(w))
    print("w0: {0}".format(w.subs({r:0})))
    Mr = -D*(w.diff(r,2) + nu/r*w.diff(r,1))
    print("Mr: {0}".format(Mr))
    Mt = -D*(w.diff(r,2) + 1/r*w.diff(r,1))
    print("Mt: {0}".format(Mt))
    Sr = 6*Mr/t**2
    print("Sr: {0}".format(Sr))
    St = 6*Mt/t**2
    print("St: {0}".format(St))
    Vmis = sympy.sqrt(Sr**2 + St**2 - Sr*St)
    print("vmis: {0}".format(Vmis))

def ex2():
    
    ''' DKT esimerkki 1, pyörähdyssymmetrinen laatta jakautuneella kuormalla '''
    
    
    t = 25e-3
    E = 100e9
    nu = 0.3
    q = 100e3
    a = 1
    D = E*t**3/(12*(1-nu**2))
    w = lambda r: -a**4*q*(1 - r**2/a**2)**2/(64*D)
    vmis = lambda r: (36*D**2*(a**2*q*(1 - r**2/a**2)/(8*D) - q*r**2/(8*D))**2/t**4 - 36*D**2*(a**2*q*(1 - r**2/a**2)/(8*D) - q*r**2/(8*D))*(a**2*nu*q*(1 - r**2/a**2)/(16*D) + a**2*q*(1 - r**2/a**2)/(16*D) - q*r**2/(8*D))/t**4 + 36*D**2*(a**2*nu*q*(1 - r**2/a**2)/(16*D) + a**2*q*(1 - r**2/a**2)/(16*D) - q*r**2/(8*D))**2/t**4)**(1/2)

    mesh = teefem.mesh.unitcircle(R=1)
    
    mat = teefem.materials.elastic(E = E, nu = nu)
    mdl = DKT(mesh = mesh)
    
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
    
    dymin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])
    
    print("DZ (acc)  : %0.14E"%(-w(0)))
    print("DZ        : %0.14E"%(dymin))
#    print("DZ (CA)   : %0.14E"%(-1.74644966293971E-01))
    
    import matplotlib.pylab as plt
#    mdl.plot()
    plotmdl(mdl, vmis)
    plt.show()

if __name__ == '__main__':
    #ex1()
    ex2()
    #rotsy()
