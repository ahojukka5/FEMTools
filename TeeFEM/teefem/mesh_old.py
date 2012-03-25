# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 17:32:59 2012

Mesh tools

Erilaisia yksinkertaisia verkotustyökaluja joilla voi muodostaa perusverkkoja
testitarkoituksiin.

@author: Jukka Aho

"""

from __future__ import division
import teefem
from teefem.common import IJV
import numpy as np
from numpy import sin,cos,pi,matrix,zeros
import matplotlib.pylab as plt
#from matplotlib import tri # Tarvii uudehkon matplotlibin
import matplotlib.delaunay
np.set_printoptions(precision = 2, edgeitems=20)


def unitsquare2(L=1,H=1):
    ''' 
    Suorakaiteen muotoinen yksinkertainen käsintehty verkotus

    P4---------GA4---------P3
    |                       |
    |                       |
    |     OM1               |
    |                       |
    |                       |
   GA1         P5          GA3
    |                       |
    |                       |
    |                       |
    |                       |
    |                       |
    P1---------GA2---------P2

    '''

    P1 = teefem.geom.Node(x=0, y=0, z=0.0)
    P2 = teefem.geom.Node(x=L, y=0, z=0.0)
    P3 = teefem.geom.Node(x=L, y=H, z=0.0)
    P4 = teefem.geom.Node(x=0, y=H, z=0.0)
    P5 = teefem.geom.Node(x=L/2, y=H/2, z=0.0)
    g1 = teefem.geom.Tria3(name = 'M1', nodes = (P1,P2,P5))
    g2 = teefem.geom.Tria3(name = 'tr = tri.Triangulation(data[:,0], data[:,1])M2', nodes = (P2,P3,P5))
    g3 = teefem.geom.Tria3(name = 'M3', nodes = (P3,P4,P5))
    g4 = teefem.geom.Tria3(name = 'M4', nodes = (P4,P1,P5))
    mesh = teefem.geom.Mesh()
    mesh.group_ma['OM1'] = set([g1,g2,g3,g4])
    mesh.group_no['GA1'] = set([P1,P4])
    mesh.group_no['GA2'] = set([P1,P2])
    mesh.group_no['GA3'] = set([P2,P3])
    mesh.group_no['GA4'] = set([P3,P4])
    mesh.group_no['P1'] = set([P1])
    mesh.group_no['P2'] = set([P2])
    mesh.group_no['P3'] = set([P3])
    mesh.group_no['P4'] = set([P4])
    mesh.group_no['P5'] = set([P5])
    P1.groups = set(['GA1','GA2','P1'])
    P2.groups = set(['GA2','GA3','P2'])
    P3.groups = set(['GA3','GA4','P3'])
    P4.groups = set(['GA4','GA1','P4'])
    g1.groups.add('OM1')
    g2.groups.add('OM1')
    g3.groups.add('OM1')
    g4.groups.add('OM1')
    mesh.nodes['N1'] = P1
    mesh.nodes['N2'] = P2
    mesh.nodes['N3'] = P3
    mesh.nodes['N4'] = P4
    mesh.nodes['N5'] = P5
    mesh.nets['M1'] = g1
    mesh.nets['M2'] = g2
    mesh.nets['M3'] = g3
    mesh.nets['M4'] = g4
    return mesh

class Trimesh(matplotlib.delaunay.Triangulation):
    ''' Kolmioverkkorutiini matplotlibin kirjastosta paranneltuna '''

    def plot(self):
        i = 0
        plt.scatter(self.x,self.y)
        for edge in self.edge_db:
            idx1 = edge[0]
            idx2 = edge[1]
            x1 = self.x[idx1]
            y1 = self.y[idx1]
            x2 = self.x[idx2]
            y2 = self.y[idx2]
            xc = 0.5*(x1+x2)
            yc = 0.5*(y1+y2)
            plt.plot([x1,x2],[y1,y2], '-b')
            plt.text(xc,yc,"%d"%(i))
            plt.text(x1,y1,idx1)
            plt.text(x2,y2,idx2)
            i += 1

    def sgm(self, fixedn):
        ''' Ei toimi. '''
        nedges = len(self.edge_db) # Reunojen määrä
        nnodes = len(self.x) # Solmupisteiden määrä
        fnodes = nnodes-fixedn # Vapaiden solmupisteiden määrä
        self.Dx = np.zeros(nedges) # "Sauvajäykkyydet"
        self.Dy = np.zeros(nedges)

        PI = 0

        lavg = 0        
        
        for edge_id in range(nedges):
            edge = self.edge_db[edge_id]
            x1 = self.x[edge[0]]
            x2 = self.x[edge[1]]
            y1 = self.y[edge[0]]
            y2 = self.y[edge[1]]
            lavg += np.sqrt((x1-x2)**2 + (y1-y2)**2)
        
        lavg /= nnodes
        print("Lagv: {0}".format(lavg))
            

        for edge_id in range(nedges):
            edge = self.edge_db[edge_id]
            x1 = self.x[edge[0]]
            x2 = self.x[edge[1]]
            y1 = self.y[edge[0]]
            y2 = self.y[edge[1]]
            l = np.sqrt((x1-x2)**2 + (y1-y2)**2)
            a = np.arctan2(y2-y1,x2-x1)
            dx = l*np.cos(a)
            dy = l*np.sin(a)
            lx = abs(self.x[edge[0]] - self.x[edge[1]])
            ly = abs(self.y[edge[0]] - self.y[edge[1]])
            D = lavg - l
            self.Dx[edge_id] = D
            self.Dy[edge_id] = D
            PI += D*l**2

        print("PI: {0}".format(PI))
        # print self.Dx
        # print self.Dy

        ijv = IJV()
        R = np.zeros(fnodes*2)
        for idx in range(fnodes):
            # Haetaan kaikki "sauvat" jotka liittyvät solmupisteeseen
            nodeid = idx+fixedn
            # print ("Nodeid: {0}:".format(nodeid))
            conlist = [nodeid in edge for edge in self.edge_db]
            # print("Connection list: {0}".format(conlist))
            #gen = (i for i,x in enumerate(conlist) if x)
            #for i in gen:
            #    print i
            edge_id_list = [i for i, x in enumerate(conlist) if x == True]
            # print("Edge id list: {0}".format(edge_id_list))
            for edge_id in edge_id_list:
                # print idx
                # print self.D[edge_id]
                ijv.add(2*idx, 2*idx, self.Dx[edge_id])
                ijv.add(2*idx+1, 2*idx+1, self.Dy[edge_id])
                edge = self.edge_db[edge_id]
                neighbor_node_id = edge[0] if edge[0] != nodeid else edge[1]
                # print("Neighbor_node_id = {0}".format(neighbor_node_id))
                if neighbor_node_id < fixedn:
                    # Tunnettuja siirtymiä
                    R[2*idx] += self.x[neighbor_node_id] * self.Dx[edge_id]
                    R[2*idx+1] += self.y[neighbor_node_id] * self.Dy[edge_id]
                    # print("Fixed")
                else:
                    # print("Free")
                    nodedim = neighbor_node_id - fixedn
                    # print("Nodedim: {0}".format(nodedim))
                    ijv.add(2*idx, 2*nodedim, -self.Dx[edge_id])
                    ijv.add(2*idx+1, 2*nodedim+1, -self.Dy[edge_id])
        #print ijv.tocoo().todense()
        #print R
        K = ijv.tocoo().tocsr()
        U = teefem.common.solve(K,R)
        print U

        #print self.x
        #print self.y

        for idx in range(fnodes):
            nodeid = idx+fixedn
            self.x[nodeid] += U[2*idx]
            self.y[nodeid] += U[2*idx+1]
    
    @property
    def average_edge_lenght(self):
        L = 0
        for edge in self.edge_db:
            x1 = self.x[edge[0]]
            x2 = self.x[edge[1]]
            y1 = self.y[edge[0]]
            y2 = self.y[edge[1]]
            L =+ np.sqrt((x1-x2)**2 + (y1-y2)**2)
        return L/len(self.edge_db)

def refine(mesh, **kwds):
    '''
    Yksinkertainen refine looppi. Katkaistaan liian pitkät reunat keskeltä.
    '''
    
    max_edge_length = kwds.get('max_edge_length', 1)
    itermax = kwds.get('iter_max', 20)
    i = 0
    breakcondition = False
    newmesh = None
    refined = True
    while not breakcondition or not refined:
        i += 1
        if i >= itermax:
            breakcondition = True
        refined = False
        for edge in mesh.edge_db:
            idx1 = edge[0]
            idx2 = edge[1]
            x1 = mesh.x[idx1]
            y1 = mesh.y[idx1]
            x2 = mesh.x[idx2]
            y2 = mesh.y[idx2]
            length = np.sqrt((x2-x1)**2 + (y2-y1)**2)
            if length > max_edge_length:
                refined = True
                mesh.x = np.append(mesh.x, 0.5*(x1+x2))
                mesh.y = np.append(mesh.y, 0.5*(y1+y2))
        newmesh = Trimesh(mesh.x,mesh.y)
    return newmesh


def unitinterval(a=0,b=1,lc=0.1):
    '''  '''
    x = [a,b]
    y = [0,0]
    x.append(0.3)
    y.append(0.6)

    nodes = len(x)
    dof = nodes*2
    lock = np.ones(2*dof)
    lock[0:4] = 0

    forcemoc(x,y,lock,plot=False,filename='output/unitinterval')


def unitsquare(a=1,b=1,lc=0.2):
    ''' Suorakaiteen muotoinen alue '''
    x = [0,a,a,0,np.random.rand(),np.random.rand(),np.random.rand()]
    y = [0,0,b,b,np.random.rand(),np.random.rand(),np.random.rand()]
        
    mesh = Trimesh(x,y)

    x0 = mesh.x
    y0 = mesh.y
#    edge_db = mesh.edge_db
    nodes = len(x0)
    dof = nodes*2
    lock = np.ones(2*dof)
    lock[0:8] = 0

    forcemoc(x,y,lock,plot=True,filename='output/unitsquare')

    return 0


def unitcircle(R=1,lc=0.5):
    ''' Yksikköympyrä '''

    x = [0]
    y = [0]

    # Ulkokehän pisteet
    outer_points = int(2*pi/lc)
    for phi in np.linspace(0,2*pi,outer_points,endpoint=False):
        x.append(R*cos(phi))
        y.append(R*sin(phi))

    # Satunnaisjakauma sisälle
    k = 1 # Kokemusperäinen vakio
    A = (lc*4**(0.25)/2)**2 # Yhden kolmion pinta-ala
    n = np.pi*R**2/A*k
    r_ = np.random.uniform(0,R,n)
    phi_ = np.random.uniform(0,2*np.pi,n)
    for (r,phi) in zip(r_,phi_):
        x.append(r*cos(phi))
        y.append(r*sin(phi))

    mesh = Trimesh(x,y)

    mesh.plot()
    plt.show()

#    mesh = refine(mesh, max_edge_length = lc)
    
    mesh.plot()
    plt.show()

    x0 = mesh.x
    y0 = mesh.y
    nodes = len(x0)
    dof = nodes*2

    lock = np.ones(2*dof)
    nn = 0
    for i in range(nodes):
        x = x0[i]
        y = y0[i]
        r = np.sqrt(x**2+y**2)
        tol = 0.1
#        tol = 0.01
        if abs(R-r) < tol:
            lock[2*i] = 0
            lock[2*i+1] = 0
            nn += 1

    print("Locked {0} nodes".format(nn))
    mesh = forcemoc(x0,y0,lock,plot=True,filename='output/unitcircle')
#    mesh.plot()
#    plt.show()
    return mesh

def sgmtest():
    ''' Testataan SGM rutiinia '''
    import sympy
    x0,x1,x2,x3,x4,y0,y1,y2,y3,y4 = sympy.var('x0,x1,x2,x3,x4,y0,y1,y2,y3,y4')
    D0,D1,D2,D3,D4,D5,D6,D7 = sympy.var('D0,D1,D2,D3,D4,D5,D6,D7')

    R0 = sympy.sqrt((x0-x1)**2 + (y0-y1)**2)
    R1 = sympy.sqrt((x1-x2)**2 + (y1-y2)**2)
    R2 = sympy.sqrt((x3-x2)**2 + (y3-y2)**2)
    R3 = sympy.sqrt((x0-x3)**2 + (y0-y3)**2)
    R4 = sympy.sqrt((x0-x4)**2 + (y0-y4)**2)
    R5 = sympy.sqrt((x1-x4)**2 + (y1-y4)**2)
    R6 = sympy.sqrt((x2-x4)**2 + (y2-y4)**2)
    R7 = sympy.sqrt((x3-x4)**2 + (y3-y4)**2)

    PI = D0*R0**2 + D1*R1**2 + D2*R2**2 + D3*R3**2 + D4*R4**2 + D5*R5**2 + D6*R6**2 + D7*R7**2
    print PI.diff(x4).expand()/2
    print PI.diff(y4).expand()/2


def force_eq(x0,y0,edge_db,lock,plot=True,filename='pomo'):
    
    ''' 
        Tasapainottaa verkon voimatasapainolla
    '''    
    
    nodes = len(x0)
    dof = nodes*2
    nedges = len(edge_db)

    # Reunan keskipituus
    Ltotal = 0
    i = 0
    for ni,nj in edge_db:
        Ltotal = np.sqrt((x0[ni]-x0[nj])**2 + (y0[ni]-y0[nj])**2)
        i += 1
    Lavg = Ltotal/nedges
    
    # Kokemus osoittaa, että stabiilissa systeemissä olisi hyvä olla vähän vetoa tai 
    # puristusta, joten säädetään sitä tällä
    k = 1.2
    Lavg *= k

    print("Lavg: {0}".format(Lavg))
    print("x0: {x0}".format(x0=x0))
    print("y0: {y0}".format(y0=y0))
    print("connections: {edge_db}".format(edge_db=edge_db))
 
    from scipy.integrate import ode
    

    def g(t,y,debug=False):
        
        def dprint(msg):
            if not debug:
                return 0
            print msg

        dprint("New iteration")
        
        n = len(y)/2
        u = y[0:n]
        ux = u[0::2]
        uy = u[1::2]

        m = np.ones(dof)*1
        dprint(np.matrix(np.diag(m)))

        K = np.matrix(np.zeros((dof,dof)))
        M = np.matrix(np.diag(m)) # + V1.T*m1*V1 + V2.T*m2*V2 
        invM = np.linalg.inv(M)

        # Vaimennusmatriisi
        alpha = 1
        beta = 0
        
        C = alpha*M + beta*K
        
        A = np.matrix(np.zeros((2*dof,2*dof)))
        A[0:dof,dof:2*dof] = np.diag(np.ones(n))
        A[dof:2*dof,0:dof] = -invM*K
        A[dof:2*dof,dof:2*dof] = -invM*C
        
        def getForceEq(nodeid):
            ''' Palauttaa voimaresultantin solmupisteessä i '''
            x = x0[nodeid] + ux[nodeid]
            y = y0[nodeid] + uy[nodeid]
            dprint("node: {node}   x: {nodex}   y: {nodey}".format(node = nodeid, nodex = x, nodey = y))
            
            # Solmupisteen yhdistävät viivalementit
            conlist = [nodeid in edge for edge in edge_db]
            edge_id_list = [i for i, iterator in enumerate(conlist) if iterator == True]
            
            Fx_tot = 0
            Fy_tot = 0
            
            for edge_idx in edge_id_list:
                edge = edge_db[edge_idx]
                neighbor_node_id = edge[0] if edge[0]!=nodeid else edge[1]
                xn = x0[neighbor_node_id] + ux[neighbor_node_id]
                yn = y0[neighbor_node_id] + uy[neighbor_node_id]
#                print("node:          {node}   x: {nodex}   y: {nodey}".format(node = nodeid, nodex = x, nodey = y))
#                print("neighbor node: {node}   x: {nodex}   y: {nodey}".format(node = neighbor_node_id, nodex = xn, nodey = yn))
                L = np.sqrt((x-xn)**2 + (y-yn)**2)
#                print("L(edge={edge}) = {length}".format(edge=edge_idx,length=L))
                dL = L-Lavg
#                print("dL(edge={edge}) = {length}".format(edge=edge_idx,length=dL))
                alp = np.arctan2(yn-y,xn-x)
                dprint("alpha(edge={edge}) = {alpha}".format(edge=edge_idx,alpha=alp*180/np.pi))
                k = 1
                F = k*dL
                Fx = F*np.cos(alp)
                Fy = F*np.sin(alp)
                dprint("F(edge={edge}) : Fx,Fy = {Fx:5.2f} , {Fy:5.2f} . Resultant: {Force}".format(edge=edge_idx,Force=F,Fx=Fx,Fy=Fy))
                Fx_tot += Fx
                Fy_tot += Fy
                
            #print("Fresultant: {0}".format(np.sqrt(Fx**2 + Fy**2)))
            
            dprint("Fx: {Fx:5.2f}  Fy: {Fy:5.2f}".format(Fx=Fx_tot,Fy=Fy_tot))
            return Fx_tot,Fy_tot
        
        # Voimavektori
        F = np.zeros(dof)
        for node_id in range(nodes):
            if lock[2*node_id]==0 or lock[2*node_id+1]==0:
                continue
            Fx,Fy = getForceEq(node_id)
            F[2*node_id] = Fx
            F[2*node_id+1] = Fy
        
        F = np.matrix(F).T
        b = np.matrix(np.zeros(2*dof)).T
        dprint("invM*F")
        dprint(invM*F)
        b[dof:2*dof,0] = invM*F
        dprint("b")
        dprint(b)

        yy = np.matrix(y).T
        dprint(yy)
        g = np.array(A*yy + b).flatten()
        dprint("g: {0}".format(g))
        dprint("lock: {0}".format(lock))
        g *= lock
        dprint("g (new): {0}".format(g))
        return g.T
    
   
    r = ode(g)
#    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(np.zeros(2*dof))
    #r.set_f_params([])
    t1 = 10
    dt = 0.1
    t = []
    y = []
    vkin = []
    condition = True
    i = 0
    itermax = 100
##    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y.real)
        v = r.y[dof:2*dof].real
        m = 1
        Vkin = np.sum(1/2*m*v**2)
        vkin.append(Vkin)
        print r.t, Vkin#, r.y.real
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False

    import matplotlib.pylab as plt
    plt.subplot(131)
    plt.plot(t,vkin,'--r')


#    u = np.array(y)
#    x1 = u[:,0] + x0[0]
#    y1 = u[:,1] + y0[0]
#    x2 = u[:,2] + x0[1]
#    y2 = u[:,3] + y0[1]
#    x3 = u[:,4] + x0[2]
#    y3 = u[:,5] + y0[2]
#    plt.subplot(132)
#    plt.plot(t,x3)
#    plt.subplot(133)
#    plt.plot(t,y3)
#    plt.show()
#    

    if plot:
        print("Plotting")
        import matplotlib.pylab as plt
        i = 0
        for u in y:
            fig = plt.figure()
    #        print("u (full)")
    #        print u
            u = u[0:dof]
    #        print("u")
    #        print u
            x = u[0::2] + x0
            y = u[1::2] + y0
            for edge in edge_db:
                x1 = x[edge[0]]
                y1 = y[edge[0]]
                x2 = x[edge[1]]
                y2 = y[edge[1]]
                plt.plot([x1,x2],[y1,y2],'--bo')
            fig.savefig('%s_%d.png'%(filename,i))
            i += 1
#        plt.show()
        print("Done")
        
#    plt.show()

#    import matplotlib.pylab as plt
#    plt.subplot(121) # X
#    plt.plot(t,x1,'--',t,x2,'--',t,x3,'--')
#    plt.subplot(122) # Y
#    plt.plot(t,y1,'--',t,y2,'--',t,y3,'--')
#    plt.show()

    return r


def forcemoc(x0,y0,lock,plot=True,filename='moc'):
    
    ''' 
        Tasapainottaa verkon voimatasapainolla.
        Tämä versio perustuu siihen että kaikki partikkelit vaikuttavat
        kaikkiin partikkeleihin. Tässä versiossa ei siis tarvitse kolmiointia
        tai muutakaan.
    '''    
    
    nodes = len(x0)
    dof = nodes*2

    print("x0: {x0}".format(x0=x0))
    print("y0: {y0}".format(y0=y0))
 
    from scipy.integrate import ode
    

    def g(t,y,debug=False):
        
        def dprint(msg):
            if not debug:
                return 0
            print msg

        dprint("New iteration")
        
        u = y[0:dof]
        ux = u[0::2]
        uy = u[1::2]

        m = np.ones(dof)*1 # Partikkeleiden massa
        dprint(np.matrix(np.diag(m)))

        K = np.matrix(np.zeros((dof,dof)))
        M = np.matrix(np.diag(m))
        invM = np.linalg.inv(M)

        # Vaimennusmatriisi
        alpha = 1
        beta = 0
        
        C = alpha*M + beta*K
        
        A = np.matrix(np.zeros((2*dof,2*dof)))
        A[0:dof,dof:2*dof] = np.diag(np.ones(dof))
        A[dof:2*dof,0:dof] = -invM*K
        A[dof:2*dof,dof:2*dof] = -invM*C
        
        def getForceEq(nodeid):
            ''' Palauttaa voimaresultantin solmupisteessä i '''
            x = x0[nodeid] + ux[nodeid]
            y = y0[nodeid] + uy[nodeid]
            dprint("node: {node}   x: {nodex}   y: {nodey}".format(node = nodeid, nodex = x, nodey = y))
            
            def dist(n1,n2):
                x1 = x0[nodeid] + ux[nodeid]
                y1 = y0[nodeid] + uy[nodeid]
                x2 = x0[nnodeid] + ux[nnodeid]
                y2 = y0[nnodeid] + uy[nnodeid]
                return np.sqrt((x1-x2)**2 + (y1-y2)**2)

            def alpha(n1,n2):
                x1 = x0[nodeid] + ux[nodeid]
                y1 = y0[nodeid] + uy[nodeid]
                x2 = x0[nnodeid] + ux[nnodeid]
                y2 = y0[nnodeid] + uy[nnodeid]
                return np.arctan2(y2-y1,x2-x1)
            
            tol = 1e-6
            
            L = np.array([dist(nodeid,nnodeid) for nnodeid in xrange(nodes)])
            Lavg = np.average(L)
            Lc = L>tol
            alp = np.array([alpha(nodeid,nnodeid) for nnodeid in xrange(nodes)])
            dL = L-Lavg
            Lc2 = dL>0
            conns = len(L)
            G = 1.0
            b = 0.5
#            for i in range(conns):
#                if abs(L[i]) < tol:
##                    print "tol: L[{0}] = {1} < {2}".format(i,L[i],tol)
#                    continue
#                val = 0#(G/L[i])**b
#                k[i] += val
            #print k

#            print("k:  {0}".format(k))
#            print("L:  {0}".format(L))
#            print("dL: {0}".format(dL))
            k = np.ones(conns)
            F = k*dL*(L>tol)*(dL<0)
#            for i in range(conns):
#                distL2 = L[i]**2
#                if distL2 < tol:
#                    continue
#                F[i] /= distL2
            Fx = F*np.cos(alp)
            Fy = F*np.sin(alp)
            Fx_tot = np.sum(Fx)
            Fy_tot = np.sum(Fy)
            
            R = 1
            r = np.sqrt(x**2 + y**2)
            if r > R:
                alp = np.arctan2(y,x)
                Fx += np.cos(alp)*(R-r)*1e6
                Fy += np.sin(alp)*(R-r)*1e6
            
#            print("L: {0}".format(L))
#            print("Lavg: {0}".format(Lavg))
#            print("Lc: {0}".format(Lc))
#            print("alp: {0}".format(alp))
#            print("dL: {0}".format(dL))
#            print("x: {x:5.2f}  y: {y:5.2f}  Fx: {Fx:5.2f}  Fy: {Fy:5.2f}  Fres: {Fres}".format(x=x,y=y,Fx=Fx_tot,Fy=Fy_tot,Fres=np.sqrt(Fx_tot**2 + Fy_tot**2)))
            return Fx_tot,Fy_tot
        
        # Voimavektori
        F = np.zeros(dof)
        for node_id in range(nodes):
            Fx,Fy = getForceEq(node_id)
            F[2*node_id] = Fx
            F[2*node_id+1] = Fy

        F = np.matrix(F).T
        b = np.matrix(np.zeros(2*dof)).T
        dprint("invM*F")
        dprint(invM*F)
        b[dof:2*dof,0] = invM*F
        dprint("b")
        dprint(b)

        yy = np.matrix(y).T
        dprint(yy)
        g = np.array(A*yy + b).flatten()
        dprint("g: {0}".format(g))
        dprint("lock: {0}".format(lock))
        g *= lock
        dprint("g (new): {0}".format(g))
        return g.T
    
   
    r = ode(g)
#    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(np.zeros(2*dof))
    #r.set_f_params([])
    t1 = 10
    dt = 0.1
    t = []
    y = []
    vkin = []
    condition = True
    i = 0
    itermax = 1000
    timemax = 300
    import time
    t0 = time.clock()
##    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y.real)
        v = r.y[dof:2*dof].real
        m = 1
        Vkin = np.sum(1/2*m*v**2)
        vkin.append(Vkin)
        print r.t, Vkin#, r.y.real
        i += 1
        if r.t >= t1:
            print("Con1")
            condition = False
        if i >= itermax:
            print("Con2")
            condition = False
        if time.clock()-t0 > timemax:
            print("Con3")
            condition = False

    import matplotlib.pylab as plt
    plt.figure(1)
    plt.plot(t,vkin,'--r')

    u = np.array(y)

    u = y[-1]
    u = u[0:dof]
    xfin = u[0::2] + x0
    yfin = u[1::2] + y0
    mesh = Trimesh(xfin,yfin)

    plt.figure(2)
    plt.scatter(xfin,yfin)
    plt.figure(3)
    mesh.plot()
    #fig.savefig('%s_%d.png'%(filename,i))

    plt.show()

    if plot:
        print("Plotting")
        i = 0
        for u in y:
            fig = plt.figure()
            ax = fig.add_subplot(111)
#            print u
#            print dof
#            print u[0:dof]
            u = u[0:dof]
            x = u[0::2] + x0
            y = u[1::2] + y0
            print x
            print y
            ax.scatter(x,y)
            fig.savefig('%s_%d.png'%(filename,i))
            i += 1
        print("Done")

    return mesh

def test1():
    x0 = np.array([np.random.rand(), 1.0, 0.0, -1.0,  0.0])
    y0 = np.array([np.random.rand(), 0.0, 1.0,  0.0, -1.0])
    edge_db = [[0,1],[0,2],[0,3],[0,4]]
    nodes = 5
    dof = nodes*2
    lock = np.ones(dof*2)
    lock[2:dof] = 0
    force_eq(x0,y0,edge_db,lock,plot=True,filename='output/test1')


if __name__ == '__main__':
    #unitsquare()
    unitcircle()
    #unitinterval()
    #sgmtest()
    #trimeshtest()
    #test1()