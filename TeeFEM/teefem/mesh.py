# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 17:32:59 2012

@author: Jukka Aho

Mesh tools

Erilaisia yksinkertaisia verkotustyökaluja joilla voi muodostaa perusverkkoja
testitarkoituksiin.

Mikään näistä ei toimi. Tästä voisi yrittää raapia kasaan jotakin toimivaa. Perusajatus
olisi että esim.

>>> import teefem
>>> mesh1 = teefem.mesh.unitcircle(r=1)
>>> print mesh1


"""

from __future__ import division
import teefem
from teefem.common import IJV
import numpy as np
from numpy import sin,cos,pi,matrix,zeros,sqrt,sum
import matplotlib.pylab as plt
#from matplotlib import tri # Tarvii uudehkon matplotlibin
import matplotlib.delaunay
np.set_printoptions(precision = 2, edgeitems=20)

def unitcircle(R=1,lc=1,meshtype='Tria3'):
    ''' 
    
    Yksikköympyrä
    
    Kolmioelementit
    
    Kentän elementit ovat elementtiryhmässä OM1 ja solmupisteryhmässä OM1
    Reuna on elementtiryhmässä GA1 ja solmupisteryhmässä GA1
    
    '''
    
    meshfile = teefem.os.path.join(teefem.datadir, 'tria3.mail') # Linear triangle mesh file
    mesh = teefem.geom.mesh(filename = meshfile)
    return mesh

###############################################################################
###############################################################################

def unitsquare(x1=0,y1=0,x2=1,y2=1,lx=10,ly=10,meshtype='Tria3'):
    ''' 
    
    Yksikkösuorakaide
    
    Ryhmät
    

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

    # Verkotusalgoritmit ovat tällä hetkellä rikki, käytetään valmista
    # verkotusdataa
    meshfile = teefem.os.path.join(teefem.datadir, 'unitsquare_tria3.msh')
    mesh = teefem.geom.mesh(filename = meshfile)
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
    
    max_edge_length = kwds.get('max_edge_length', 0.1)
    itermax = kwds.get('iter_max', 100)
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


def unitsquare2(x1,y1,x2,y2,lx,ly,meshtype='Tria3'):
    ''' 
        Suorakaiteen muotoinen alue
    '''
    x = np.linspace(x1,x2,lx)
    y = np.linspace(y1,y2,ly)
    p = np.array([[xi,yi] for xi in x for yi in y])
    X = p[:,0]
    Y = p[:,1]
    print X
    print Y
    
    if meshtype != 'Tria3':
        return None
    trmesh = Trimesh(X,Y)
    trmesh.plot()
    plt.show()
    print("Kesken")
    return None
    
    mesh = teefem.geom.Mesh()

#    mesh.nodes = set()
#    for xi,yi in zip(trmesh.x, trmesh.y):
#        node = teefem.geom.Node(x = xi, y = yi)
#        mesh.nodes.add(node)
    
    

    return mesh
#    mesh = refine(mesh, max_edge_length = 0.1)
#    
#    mesh.plot()
#    plt.show()

#    x0 = mesh.x
#    y0 = mesh.y
#    edge_db = mesh.edge_db
#
#    nodes = len(x0)
#    dof = 2*nodes
#    dd = len(xa)*4
#    lock = np.ones(2*dof)
#    tol = 1e-5
#    # 0 jos lukittu, 1 jos vapaa
#    for i in range(nodes):
#        if abs(x0[i]-a) < tol or abs(x0[i]+a) < tol:
#            lock[2*i] = 0
#        if abs(y0[i]-b) < tol or abs(y0[i]+b) < tol:
#            lock[2*i+1] = 0
#    lock[0:dd] = 0
#
#    force_eq(x0,y0,edge_db,lock,plot=True,filename='output/unitsquare')

    return 0


def unitcircle1(R=1,lc=0.5):
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

#    mesh.plot()
#    plt.show()

#    mesh = refine(mesh, max_edge_length = lc)
#    mesh.plot()
#    plt.show()

    x0 = mesh.x
    y0 = mesh.y
    
    distance_function = lambda x,y: x**2+y**2-1

    mesh = force_eq(x0,y0,mesh.edge_db,distance_function,plot=True,filename='output/unitcircle')
#    mesh.plot()
#    plt.show()
    return mesh

def unitcircle2(R=1,lc=1.0):
    ''' Toinen versio yksikköympyrästä, alkujakauma kuten distmeshissä '''
    
    fd = lambda x,y: x**2 + y**2 - R
    fh = lambda x,y: 1
    pfix = [(0,0)] # Halutaan keskipiste

    # 1. Create initial distribution in bounding box (equilateral triangles)
    xx = np.arange(-R,R,lc)
    yy = np.arange(-R,R,lc*sqrt(3)/2)
    x,y = np.meshgrid(xx,yy)
    x[2::,:] += 0.5*lc # Shift even rows
    px = x.flatten()
    py = y.flatten()
    p = np.array(zip(px,py))
    
    # 2. Remove points outside the region, apply the rejection method
    p = np.array([pi for pi in p if fd(*pi) < 0])
    r0 = 1/np.array([fh(*pi) for pi in p])**2 # Probability to keep point
    r0 /= max(r0) # Normeerataan
    p = np.array([ p[i] for i in xrange(len(r0)) if np.random.rand() < r0[i] ])
    p = np.concatenate((p,pfix),axis=0)
    
    x0 = [pi[0] for pi in p]
    y0 = [pi[1] for pi in p]
    mesh = Trimesh(x0,y0)
    
    return force_eq(mesh,fd,fh)

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


def force_eq(mesh,fd,fh,plot=True,filename='output/force_eq'):
    
    ''' 
        Tasapainottaa pistejoukon voimatasapainolla käyttäen kolmiointia.
        
    '''

    def g(t,y,args):
        
        mesh = args[0]
        edge_db = mesh.edge_db
        x0 = mesh.x
        y0 = mesh.y
        fd = args[1]
        fh = args[2]
        
        n = len(y)/2
        u = y[0:n]
        ux = u[0::2]
        uy = u[1::2]

        m = np.ones(dof)*1

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
        
        def getForceEq():
            ''' Palauttaa voimaresultantin '''

            def dist(id1,id2):
                x1 = x0[id1] + ux[id1]
                y1 = y0[id1] + uy[id1]
                x2 = x0[id2] + ux[id2]
                y2 = y0[id2] + uy[id2]
                return np.sqrt((x2-x1)**2+(y2-y1)**2)
            
            def angle(id1,id2):
                x1 = x0[id1] + ux[id1]
                y1 = y0[id1] + uy[id1]
                x2 = x0[id2] + ux[id2]
                y2 = y0[id2] + uy[id2]
                return np.arctan2(y2-y1,x2-x1)

            def mid(id1,id2):
                x1 = x0[id1] + ux[id1]
                y1 = y0[id1] + uy[id1]
                x2 = x0[id2] + ux[id2]
                y2 = y0[id2] + uy[id2]
                return (0.5*(x1+x2), 0.5*(y1+y2))

            Fscale = 1.2
            # Sauvojen pituudet
            L = np.array([dist(id1,id2) for (id1,id2) in edge_db])
            
            # Sauvojen keskipisteet
            midpoints = np.array([mid(id1,id2) for (id1,id2) in edge_db])
            
            # Sauvojen "h-funktiot"
            hbars = np.array([fh(*xy) for xy in midpoints])
            
            # Sauvojen haluttu pituus
            L0 = hbars * Fscale * np.sum(L)/np.sum(hbars)
            
            # Sauvavoimat 2.5
            k = 1
            F = k*(L0-L)/L
            F *= (F>0) # Vain puristava voima: kerrotaan ykkösellä kaikki positiiviset voimat, negatiiviset nollalla
            
#            print("L: {0}".format(L))
#            print("midpoints: {0}".format(midpoints))
#            print("hbars: {0}".format(hbars))
#            print("L0: {0}".format(L0))
#            print("F: {0}".format(F))

            Fvec = np.zeros(dof)

            nedges = len(edge_db)

            for edge_idx in xrange(nedges):
                
                edge = edge_db[edge_idx]
                
                id1 = edge[0]
                id2 = edge[1]
                
                x1 = x0[id1] + ux[id1]
                y1 = y0[id1] + uy[id1]
                x2 = x0[id2] + ux[id2]
                y2 = y0[id2] + uy[id2]

                l = L[edge_idx]
                Fi = F[edge_idx]
                dy = y2-y1
                dx = x2-x1

                cosa = dx/l
                sina = dy/l
                
                u1 = 2*id1
                u2 = 2*id1+1
                u3 = 2*id2
                u4 = 2*id2+1
                
                Fix = Fi*cosa
                Fiy = Fi*sina

#                print("Bar id {id}: Fi = {Fix:5.2f}i + {Fiy:5.2f}j = |{Fi:5.2f}|".format(id = edge_idx, Fi = Fi, Fix = Fix, Fiy=Fiy))
                
                Fvec[u1] -= Fix
                Fvec[u2] -= Fiy
                Fvec[u3] += Fix
                Fvec[u4] += Fiy

            FF = 1e1 # "Reunavoima"
            for node_id in xrange(nnodes):
                nodex = x0[node_id] + ux[node_id]
                nodey = y0[node_id] + uy[node_id]
                fval = fd(nodex,nodey)
                if fval > 0:
#                    print("Node {nodeid} it outside -> fd({x:5.2f},{y:5.2f}) = {fd:5.2f} > 0".format(nodeid=node_id,x=nodex,y=nodey,fd=fval))
                    dx = nodex - 0
                    dy = nodey - 0
                    l = np.sqrt(nodex**2 + nodey**2)
                    cosa = nodex/l
                    sina = nodey/l
                    Fi = fval*FF
                    Fx = -Fi*cosa
                    Fy = -Fi*sina
#                    print("Node {nodeid}: adding 'edge' force {Fx:5.2f}i + {Fy:5.2}j = {Fi:5.2}".format(nodeid = node_id, Fx = Fx, Fy = Fy, Fi = Fi))
                    Fvec[2*node_id] += Fx
                    Fvec[2*node_id+1] += Fy
            return Fvec
        
        # Voimavektori
        F = getForceEq()
        F = np.matrix(F).T
        b = np.matrix(np.zeros(2*dof)).T
        b[dof:2*dof,0] = invM*F

        yy = np.matrix(y).T
        g = np.array(A*yy + b).flatten()
        return g.T

###############################################################################
    
    from scipy.integrate import ode

    nnodes = len(mesh.x)
    dof = nnodes*2

    print("x0: {x0}".format(x0=mesh.x))
    print("y0: {y0}".format(y0=mesh.y))

    print("connections: {edge_db}".format(edge_db=mesh.edge_db))
 

    r = ode(g)
#    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(np.zeros(2*dof))
    r.set_f_params([mesh,fd,fh])

    dt = 0.1
    t = []
    vkinlist = []
    condition = True
    i = 0
    
    # Pysäytysehdot
    itermax = 1e5
    timemax = 300
    vkinmin = 1e-3 # Kinemaattinsen energian minimi simuloinnin pysäyttämiseen
    tol = 0.03 # Kinemaattisen energian maksimi uudelleenkolmiointiin 
    
    phi = np.linspace(0,2*pi,500)
    radius = 1
    cx = radius*cos(phi)
    cy = radius*sin(phi)

    import time
    time0 = time.clock()

    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)

        u = r.y[0:dof].real
        v = r.y[dof:2*dof].real

        ux = u[0::2]
        uy = u[1::2]

        x = mesh.x + ux
        y = mesh.y + uy

        m = np.ones(len(v))
        M = np.sum(m)
        vkin = np.sum(1/2*m*v**2)/M
        # Kokeillaan että aika-askel on verrannollinen liike-energiaan
        dt = max(vkin,0.03)
        if vkin>tol:
            print("vkin = {vkin:5.2f} > {tol:5.2f}. Retriangulating.".format(vkin=vkin,tol=tol))
            mesh = Trimesh(x,y)
            r.set_f_params([mesh,fd,fh])
        umax = max(np.sqrt(ux**2 + uy**2))
        utol = 100e6
        if umax>utol:
            print("umax = {umax:5.2f} > {utol:5.2f}. Retriangulating.".format(umax=umax,utol=utol))
            mesh = Trimesh(x,y)
            r.set_f_params([mesh,fd,fh])
        vkinlist.append(vkin)

        
        time1 = time.clock()
        deltat = time1-time0
        print("Real time: {rtime:5.2f} Simulation time: {stime:5.2f}  Kinetic energy (normed): {vkin:5.2f} umax: {umax:5.2f}".format(rtime = deltat, stime=r.t, vkin=vkin, umax=umax))

        if plot:
            fn = '%s_%d.png'%(filename,i)
            print("Plotting to %s"%(fn))
            fig = plt.figure(1)
            fig.clf()
            plt.plot(cx,cy,'--r')
            for edge in mesh.edge_db:
                x1 = x[edge[0]]
                y1 = y[edge[0]]
                x2 = x[edge[1]]
                y2 = y[edge[1]]
                plt.plot([x1,x2],[y1,y2],'--bo')
            plt.xlim((-1,1))
            plt.ylim((-1,1))
            fig.savefig(fn)

        i += 1

        if i >= itermax or deltat > timemax or vkin < vkinmin: 
            condition = False

    plt.figure(2)
    plt.subplot(111)
    plt.plot(t,vkinlist,'--r')
    plt.show()

#    u = np.array(y)
#    x1 = u[:,0] + x0[0]
#    y1 = u[:,1] + y0[0]
#    x2 = u[:,2] + x0[1]
#    y2 = u[:,3] + y0[1]
#    x3 = u[:,4] + x0[2]
#    y3 = u[:,5] + y0[2]
#    print x1
#    print y1
#    print x2
#    print y2
#    print x2
#    print y3
#    plt.figure(3)
#    plt.subplot(231)
#    plt.plot(t,x1)
#    plt.subplot(232)
#    plt.plot(t,y1)
#    plt.subplot(233)
#    plt.plot(t,x2)
#    plt.subplot(234)
#    plt.plot(t,y2)
#    plt.subplot(235)
#    plt.plot(t,x3)
#    plt.subplot(236)
#    plt.plot(t,y3)
#    plt.show()

    return mesh


def forcemoc(x0,y0,plot=True,filename='moc'):
    
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
                return np.arctan2(y1-y2,x1-x2)
            
            tol = 1e-6
            
            L = np.array([dist(nodeid,nnodeid) for nnodeid in xrange(nodes)])
            Fscale = 1.2
            Lavg = np.sum(L)/(len(L)-1)
            L0 = Lavg*Fscale
            dL = L0-L

            alp = np.array([alpha(nodeid,nnodeid) for nnodeid in xrange(nodes)])

#            print("k:  {0}".format(k))
#            print("L:  {0}".format(L))
#            print("dL: {0}".format(dL))
            k = 1
            F = k*dL*(L>tol)*(dL>0)
            Fx = F*np.cos(alp)
            Fy = F*np.sin(alp)
            Fx_tot = np.sum(Fx)
            Fy_tot = np.sum(Fy)
            
            R = 1
            r = np.sqrt(x**2 + y**2)
            if r > R:
                alp = np.arctan2(y,x)
                Fx_tot += np.cos(alp)*(R-r)*1e2
                Fy_tot += np.sin(alp)*(R-r)*1e2
            
#            print("L: {0}".format(L))
#            print("Lavg: {0}".format(Lavg))
#            print("alp: {0}".format(alp))
#            print("dL: {0}".format(dL))
            #print("x: {x:5.2f}  y: {y:5.2f}  Fx: {Fx:5.2f}  Fy: {Fy:5.2f}  Fres: {Fres}".format(x=x,y=y,Fx=Fx_tot,Fy=Fy_tot,Fres=np.sqrt(Fx_tot**2 + Fy_tot**2)))
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
#        dprint("g: {0}".format(g))
#        dprint("lock: {0}".format(lock))
#        g *= lock
#        dprint("g (new): {0}".format(g))
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
    itermax = 10000000
    timemax = 15*60
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

    plt.figure(2)
    
    phi = np.linspace(0,2*pi,500)
    R = 1
    plt.plot(R*np.cos(phi),R*np.sin(phi),'--b')
    
    
    plt.scatter(xfin,yfin)
#    plt.figure(3)
    #fig.savefig('%s_%d.png'%(filename,i))


    if plot:
        print("Plotting")
        i = 0
        for u in y:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            phi = np.linspace(0,2*pi,500)
            R = 1
            plt.plot(R*np.cos(phi),R*np.sin(phi),'--b')
            u = u[0:dof]
            x = u[0::2] + x0
            y = u[1::2] + y0
            print x
            print y
            ax.scatter(x,y)
            fig.savefig('%s_%d.png'%(filename,i))
            i += 1
        print("Done")


    return Trimesh(xfin,yfin)

def test1():
    x0 = np.array([np.random.rand(), 1.0, 0.0, -1.0,  0.0])
    y0 = np.array([np.random.rand(), 0.0, 1.0,  0.0, -1.0])
    edge_db = [[0,1],[0,2],[0,3],[0,4]]
    nodes = 5
    dof = nodes*2
    lock = np.ones(dof*2)
    lock[2:dof] = 0
    force_eq(x0,y0,edge_db,lock,plot=True,filename='output/test1')

def forcemoc_test1():
#    r = np.random.rand(100)
#    phi = np.random.rand(100)*2*np.pi
#    x = r*np.cos(phi)
#    y = r*np.sin(phi)
    mesh = forcemoc(x,y,plot=True,filename='output/unitcircle_moc')

def force_eq_test1():
#    x = np.array([-0.31,0.44,-0.13])
#    y = np.array([-0.53,-0.24,0.46])
    r = np.random.rand(50)
    phi = np.random.rand(50)*2*np.pi
    x = r*cos(phi)
    y = r*sin(phi)
#    m = Trimesh(x,y)
    fd = lambda x,y: x**2+y**2-1
    mesh = force_eq(x,y,distance_function = fd)

if __name__ == '__main__':
    #forcemoc_test1()
    #force_eq_test1()
    unitsquare(0,0,2,2,10,10)
    #unitcircle2()
    #unitcircle()
    #unitinterval()
    #sgmtest()
    #trimeshtest()
    #test1()