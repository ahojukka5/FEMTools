# -*- coding: utf-8 -*-
"""
Created on Sat Mar 24 02:27:29 2012

@author: Jukka

Kaikenlaisia verkotukseen liittyviä kokeiluja

"""

def trimeshtest():
    ''' Kolmioverkon testausta '''
    x = [1,0,-1,0,0.1,-0.2]
    y = [0,1,0,-1,0.1,-0.3]
    mesh = Trimesh(x,y)
    # Haetaan "voimatasapaino", kiinnitetään listan neljä ensimmäistä
    # solmupistettä
    for i in range(3):
        mesh.sgm(4)
    mesh.plot()
    plt.show()
    
def onedmeshtest():
    ''' 
    Yksidimensioinen verkon tasapainotustesti 

    Perustuu suuriin siirtymiin ja staattiseen analyysiin
    
    '''
    # Kiinnitetyt vapausasteet
    x0 = 0
    x1 = 1
    # Optimoitava parametri
    x2 = 0.4
    
    itermax = 10
    i = 0
    x4 = 0.1
    while i < itermax:
        i += 1
        # Jäykkyysmatriisit
        k1 = E*A/L1

        d4l = np.sqrt((x4-x0)**2)
        d6l = np.sqrt((x4-x2)**2)
        avgl = 0.5*(d4l+d6l)
        D4 = d4l/avgl
        D6 = d6l/avgl
        x4 += (D4*x0+D6*x2)
        print("d4l: {0:5.2f}   d6l: {1:5.2f}   avgl: {2:5.2f}   D4: {3:5.2f}   D6: {4:5.2f}   x4: {5:5.2f}".format(d4l,d6l,avgl,D4,D6,x4))


def onedmeshtest2():
    ''' Yksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = F
        
        x0|---x2------|x1
        
    '''    
    
    # Alkutila
    x0 = 0
    x1 = 1
    x2 = np.random.rand()
    Lavg = 1
    nnodes = 3 # Solmupisteiden määrä
    fixednodes = 2 # Kiinnitettyjen solmupisteiden määrä
    freenodes = nnodes - fixednodes
    nedges = 2 # Reunojen määrä
    
    edge_db = [[0,2],[1,2]]
    
    # Alkunopeudet
    v0 = np.zeros(freenodes)

    # Massa
    m = 1

    alpha = 2
    beta = 0
    
    from scipy.integrate import ode

    y0 = [x2, v0]
    t0 = 0

    E = 1
    A = 1
    k1 = E*A*Lavg
    k2 = E*A*Lavg
    M = m*np.matrix(np.ones(freenodes)) 
    invM = np.linalg.inv(M)
    #K = np.matrix([k1+k2])
    K = np.matrix([0])
    C = alpha*M + beta*K
    
    n = 1
    A = np.matrix(np.zeros((2*n,2*n)))
    A[0:n,n:2*n] = np.diag(np.ones(n))
    A[n:2*n,0:n] = -invM*K
    A[n:2*n,n:2*n] = -invM*C
    print A
    
    print("invM")
    print invM
    
    def arctan2(y,x):
        if y == 0:
            if x >= 0:
                return 0
            else:
                return np.pi
        alp = np.arctan(y/x)
        if x < 0: 
            alp = alp*(-1)
        return float(alp)
    
    def dist(x1,x2):
        return np.sqrt((x1-x2)**2)

    def const_f(y):
        # b-matriisi
        f = np.zeros(n)
        for idx in range(n):
#            print("R{0}:".format(idx))
            x2 = y[0]

            L1 = dist(x2,x0)
            dL1 = L1-Lavg
            # Fortran-tulkki ei ymmärrä np.arctan2
#            alp = np.arctan2(0,x0-x2)
            alp = arctan2(0,x0-x2)
            F = k1*dL1
            Fx = F*np.cos(alp)
            f[idx] += Fx

            L2 = dist(x2,x1)
            dL2 = L2-Lavg
            alp = arctan2(0,x1-x2)
            F = k2*dL2
            Fx = F*np.cos(alp)
            f[idx] += Fx

#            print("L1: {0:5.2f} L2: {1:5.2f} dL1: {2:5.2f} dL2: {3:5.2f} fx: {4:5.2f}".format(L1,L2,dL1,dL2,f[idx]))
            #print f[idx]
        return f
        


    def g(t,y):
        bb = invM*np.matrix(const_f(y))
        b = np.zeros(2*n)
        b[n:2*n] = bb
        b = np.matrix(b).T
        # print("b: {0}".format(b))
        g = A*np.matrix(y).T + b
        #print g
        return g
    
#    dt = 0.1
#    y0 = np.matrix([x2,0]).T
#    g0 = g(t0,y0)
#    print g0
#    g1 = g(t0+dt,g0)
#    print g1
#   
    r = ode(g)
    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(y0,t0)
#    r.set_f_params([m,c,k])
    t1 = 20
    dt = 0.1
    t = []
    y = []
    condition = True
    i = 0
    itermax = 200
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y[0].real)
        print r.t, r.y[0].real
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False
    import matplotlib.pylab as plt
    plt.plot(t,y,'--')
    plt.show()
    return r

def arctan2(y,x):
    if y == 0:
        if x >= 0:
            return 0
        else:
            return np.pi
    alp = np.arctan(y/x)
    if x < 0: 
        alp = alp*(-1)
    return float(alp)

def dist(x1,x2):
    return np.sqrt((x1-x2)**2)


def onedmeshtest3():
    ''' Yksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = F
        
        x0|---x2------|x1
        
    '''    
    
    # Alkutila
    x0 = np.array([0,1,np.random.rand()])
    edge_db = [[0,2],[2,1]]
    nnodes = len(x0) # Solmupisteiden määrä
    fixednodes = 2 # Kiinnitettyjen solmupisteiden määrä
    freenodes = nnodes - fixednodes
    nedges = len(edge_db) # Reunojen määrä
    # Keskipituus
    Ltotal = 0
    for ni,nj in edge_db:
        Ltotal += np.sqrt((x0[ni]-x0[nj])**2)
    Lavg = Ltotal/nedges
    print("Lavg: {0}".format(Lavg))
    
    # Alkunopeudet
    v0 = np.zeros(freenodes)

    # Massa
    m = 1

    alpha = 2
    beta = 0
    
    from scipy.integrate import ode

    y0 = (x0[2],0)
    t0 = 0

    E = 1
    A = 1

    k = E*A*Lavg*np.ones(nedges)


    M = m*np.matrix(np.diag(m*np.ones(freenodes)))

    invM = np.linalg.inv(M)
    #K = np.matrix([k1+k2])
    K = np.matrix(np.zeros((freenodes,freenodes)))
    C = alpha*M + beta*K
    
    n = 1
    A = np.matrix(np.zeros((2*n,2*n)))
    A[0:n,n:2*n] = np.diag(np.ones(n))
    A[n:2*n,0:n] = -invM*K
    A[n:2*n,n:2*n] = -invM*C
    print A
    
    print("invM")
    print invM
    
    
    def dist(x1,x2):
        return np.sqrt((x1-x2)**2)

    def const_f(y):
        # b-matriisi
        f = np.zeros(n)
        for idx in range(n):
            print("R{0}:".format(idx))
            nodeid = idx+fixednodes
            nodecoord = y[idx]
            for edge in edge_db:
#                print("Reuna: {0}".format(edge))
                neighbor_node_id = edge[0] if edge[0]!=nodeid else edge[1]
#                print("Naapuri: {0}".format(neighbor_node_id))
                nncoord = x0[neighbor_node_id] if neighbor_node_id < fixednodes else y[neighbor_node_id]
#                print("x: {0:5.2f} x_naapuri {1:5.2f}".format(nodecoord, nncoord))
                length = dist(nodecoord, nncoord)
                dLength = length - Lavg
                # Fortran-tulkki ei ymmärrä np.arctan2
    #            alp = np.arctan2(0,x0-x2)
                alp = arctan2(0,nncoord-nodecoord)
#                print("Suunta: {0}".format(alp*180.0/pi))
                F = k[idx]*dLength
                Fx = F*np.cos(alp)
#                print("Voima: {0}".format(Fx))
                f[idx] += Fx
#            print("Kokonaisvoima: {0}".format(f[idx]))
#            print("L1: {0:5.2f} L2: {1:5.2f} dL1: {2:5.2f} dL2: {3:5.2f} fx: {4:5.2f}".format(L1,L2,dL1,dL2,f[idx]))
            #print f[idx]
        return f
        


    def g(t,y):
        bb = invM*np.matrix(const_f(y))
        b = np.zeros(2*n)
        b[n:2*n] = bb
        b = np.matrix(b).T
        # print("b: {0}".format(b))
#        print A
#        print np.matrix(y)
#        print b
        g = A*np.matrix(y).T + b
        #print g
        return g
    
#    dt = 0.1
#    y0 = np.matrix([x2,0]).T
#    g0 = g(t0,y0)
#    print g0
#    g1 = g(t0+dt,g0)
#    print g1
#   
    r = ode(g)
    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(y0)
#    r.set_f_params([m,c,k])
    t1 = 10
    dt = 0.1
    t = []
    y = []
    condition = True
    i = 0
    itermax = 200
#    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y[0].real)
        print r.t, r.y[0].real
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False
    import matplotlib.pylab as plt
    plt.plot(t,y,'--')
    plt.show()
    return r

def onedmeshtest4():
    ''' Yksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = F
        
        x0|---x2------|x1
        
    '''    
    
    # Alkutila
    x0 = np.array([0,1,np.random.rand(),np.random.rand()])
    edge_db = [[0,2],[1,2]]
    nnodes = len(x0) # Solmupisteiden määrä
    fixed_nodes = range(2) # Kiinnitetyt solmupisteet
#    freenodes = nnodes - len(fixednodes)
    nedges = len(edge_db) # Reunojen määrä
    # Keskipituus
    Ltotal = 0
    for ni,nj in edge_db:
        Ltotal += np.sqrt((x0[ni]-x0[nj])**2)
    Lavg = Ltotal/nedges
    print("Lavg: {0}".format(Lavg))
    
    # Alkunopeudet
    v0 = np.zeros(nnodes)

    # Massa
    m = 1

    alpha = 2
    beta = 0
    
    from scipy.integrate import ode

    y0 = np.zeros(2*nnodes)
    y0[0:nnodes] = v0
    y0[nnodes:2*nnodes] = x0
    print("Y0")
    print y0

    E = 1
    A = 1

    k = E*A*Lavg*np.ones(nedges)


    M = m*np.diag(m*np.ones(nnodes))
#    M = m*np.diag((0,0,1))
    invM = np.linalg.inv(M)

    print("invM")
    print invM

    #K = np.matrix([k1+k2])
    K = np.zeros((nnodes,nnodes))
    C = alpha*M + beta*K
    print("C")
    print C
    
    n = nnodes
    A = matrix(zeros((2*n,2*n)))
    A[0:n,n:2*n] = np.diag(np.ones(n))
    A[0:n,n:2*n] = np.diag((1,1,0,0))
    A[n:2*n,0:n] = -invM*K
    A[n:2*n,n:2*n] = -invM*C
    print("A")
    print A
    
    def const_f(res):
        # b-matriisi
        f = np.zeros(n)
        for nodeid in range(n):
#            print("R{0}:".format(nodeid))
            x = res[nodeid]
            for edge_idx in range(len(edge_db)):
                edge = edge_db[edge_idx]
#                print("Reuna: {0}".format(edge))
                neighbor_node_id = edge[0] if edge[0]!=nodeid else edge[1]
#                print("Naapuri: {0}".format(neighbor_node_id))
                xn = res[neighbor_node_id]
                L = dist(x, xn)
                dL = L - Lavg
                alp = arctan2(0,xn-x)
                F = k[edge_idx]*dL
                Fx = F*np.cos(alp)
                f[nodeid] += Fx
#            if nodeid in fixed_nodes:
#                f[nodeid] = 0
#            print("Kokonaisvoima: {0}".format(f[nodeid]))
        return f
        
    def g(t,res):
        f = const_f(res)
        print("f")
        print f
        bb = np.array(np.matrix(invM)*np.matrix(f).T).flatten()
        print("bb")
        print bb
        b = np.zeros(2*n)
        b[n:2*n] = bb
        print("b")
        print b
        b = np.matrix(b).T
        y = np.matrix(res).T
#        print("b,y")
#        print b
#        print y
        g = A*y + b
        return g
    
    r = ode(g)
    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(y0)
#    r.set_f_params([m,c,k])
    t1 = 200
    dt = 0.01
    t = []
    y = []
    condition = True
    i = 0
    itermax = 2000
#    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y)
        print "{0:5.2f} : {1}".format(r.t, r.y.real)
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False
    import matplotlib.pylab as plt
    plt.plot(t,y,'--')
    plt.show()
    return r

def onedmeshtest5():
    ''' Yksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = F
        
        x0|------x1------|x2
        
    '''    
    
    # Alkutila
    x0 = np.array([0.0, 0.4, 1.0])
    y0 = np.array([0.0, 0.1, 0.0])
    edge_db = [[0,1],[1,2]]

    # Täytyy olla jokin käsitys alkutilasta. Tässä oletetaan että
    # 1) viivaelementit ovat yhtä pitkiä (keskipituus)
    # 2) venymä on poikkeama tästä keskipituudesta

    nedges = len(edge_db)
    u0 = np.zeros(nedges)
    
    # Reunan keskipituus
    Ltotal = 0
    i = 0
    for ni,nj in edge_db:
        u0[i] = np.sqrt((x0[ni]-x0[nj])**2 + (y0[ni]-y0[nj])**2)
        Ltotal += u0[i]
        i += 1
    Lavg = Ltotal/nedges

    print("Lavg: {0}".format(Lavg))
    print("x0: {x0}".format(x0=x0))
    print("y0: {y0}".format(y0=y0))
    print("u0: {u0}".format(u0=u0))
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
        v = y[n:2*n]
        dprint("y")
        dprint(y)
        dprint("u")
        dprint(u)
        dprint("v")
        dprint(v)
        dprint("ux")
        dprint(ux)
        dprint("uy")
        dprint(uy)

        # Jäykkyysmatriisi
        def k(edge = None):

            node1 = edge[0]
            node2 = edge[1]

            dprint("node1: {0}    node2: {1}".format(node1,node2))

            x1 = float(x0[node1]) + float(ux[node1])
            y1 = float(y0[node1]) + float(uy[node1])
            x2 = float(x0[node2]) + float(ux[node2])
            y2 = float(y0[node2]) + float(ux[node2])
            
            ux1 = float(ux[node1])
            uy1 = float(uy[node1])
            ux2 = float(ux[node2])
            uy2 = float(ux[node2])
            dprint("uy: {0}".format(uy))
            
            dprint("ux1: {ux1:5.2f}  uy1: {uy1:5.2f}  ux2: {ux2:5.2f}  uy2: {uy2:5.2f}".format(ux1=ux1,uy1=uy1,ux2=ux2,uy2=uy2))
            
            alp = np.arctan2(y2-y1,x2-x1)
            E = 1e6
            A = 1e-6
            EA = 50e2
            L = Lavg#np.sqrt((x2-x1)**2 + (y1-y2)**2)
            
            # "Stress stiffening" kerroin on verrannollinen venymään
            du = np.sqrt((ux2-ux1)**2 + (uy2-uy1)**2)
            dprint("du: {0} , alp: {1} , x1: {x1:5.2f} y1: {y1:5.2f} x2: {x2:5.2f} y2: {y2:5.2f}".format(du,alp,x1=x1,y1=y1,x2=x2,y2=y2))
            
            c = np.cos(alp)
            s = np.sin(alp)
            cc = c*c
            ss = s*s
            cs = c*s
            ki = EA/L*np.matrix([[cc,cs,-cc,-cs],[cs,ss,-cs,-ss],[-cc,-cs,cc,cs],[-cs,-ss,cs,ss]])

#            du = 0
            dprint("ki:")
            dprint(ki*(1+du))
            return ki*(1+du)
        
        def m(edge = None, rho = 1):
            ''' Massamatriisin rakennus '''
            x1 = float(x0[edge[0]])
            y1 = float(y0[edge[0]])
            x2 = float(x0[edge[1]])
            y2 = float(y0[edge[1]])
            L = np.sqrt((x2-x1)**2 + (y1-y2)**2)
            M = rho*L/6*np.matrix([[2,1],[1,2]])
            alp = np.arctan2(y2-y1,x2-x1)
            c = np.cos(alp)
            s = np.sin(alp)
            B = np.matrix([[c,s,0,0],[0,0,c,s]])
            return B.T*M*B

        k1 = k(edge = edge_db[0])
        k2 = k(edge = edge_db[1])
        m1 = m(edge = edge_db[0])
        m2 = m(edge = edge_db[1])
        
        V1 = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0]])
        V2 = np.matrix([[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])

        dprint(V1)
        dprint(k1)
        dprint(V2)
        dprint(k2)

        dof = 6 # Vapausasteiden määrä

        m = np.ones(dof)*1
        dprint(np.matrix(np.diag(m)))

        K = V1.T*k1*V1 + V2.T*k2*V2
        dprint("K")
        dprint(K)
        M = np.matrix(np.diag(m)) # + V1.T*m1*V1 + V2.T*m2*V2 
        dprint("M")
        dprint(M)
        
        # Sakkomatriisit
        # Näitä ei tarvita jos "nollataan" vapausasteita g-vektorista
#        a1 = matrix([1,0,0,0,0,0])
#        a2 = matrix([0,1,0,0,0,0])
#        a3 = matrix([0,0,0,0,1,0])
#        a4 = matrix([0,0,0,0,0,1])
#        s = 1e1
#        K += s*(a1.T*a1 + a2.T*a2 + a3.T*a3 + a4.T*a4)

#        # Massamatriisi
#        M = m*np.matrix(np.diag(m*np.ones(dof)))

        invM = np.linalg.inv(M)
        # Vaimennusmatriisi
        alpha = 10
        beta = 0
        
        C = alpha*M + beta*K
        
        A = np.matrix(np.zeros((2*dof,2*dof)))
        A[0:dof,dof:2*dof] = np.diag(np.ones(n))
        A[dof:2*dof,0:dof] = -invM*K
        A[dof:2*dof,dof:2*dof] = -invM*C
        dprint("A")
        dprint(A)
        
        F = np.matrix([0,0,0,-10,0,0]).T
        b = np.matrix(np.zeros(2*dof)).T
        dprint("invM*F")
        dprint(invM*F)
        b[dof:2*dof,0] = invM*F
        dprint("b")
        dprint(b)

        yy = np.matrix(y).T
        dprint(yy)
        g = np.array(A*yy + b).flatten()
        # Eliminoidaan vapausasteita vaikkapa tällä tavalla
        lock = np.ones(dof*2)
        lock[0] = 0
        lock[1] = 0
        lock[4] = 0
        lock[5] = 0
        dprint("g: {0}".format(g))
        dprint("lock: {0}".format(lock))
        g *= lock
        dprint("g (new): {0}".format(g))
        return g.T
    
#    dt = 0.1
#    t0 = 0
    dof = 6    
    y0 = np.zeros(2*dof)
    #y0[0:dof] = u0
##
#    g0 = g(t0,y0)
#    print("g0")
#    print g0
#    g1 = g(t0+dt,g0)
#    print("g1")
#    print g1
   
    r = ode(g)
#    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(y0)
##    r.set_f_params([m,c,k])
    t1 = 1
    dt = 0.01
    t = []
    y = []
    condition = True
    i = 0
    itermax = 1000
##    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y.real)
        print r.t, r.y.real
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False

    n = len(y)/2
    y = np.array(y)
    print y
    u1 = y[:,0]
    u2 = y[:,1]
    u3 = y[:,2]
    u4 = y[:,3]
    u5 = y[:,4]
    u6 = y[:,5]

    import matplotlib.pylab as plt
    plt.subplot(121) # Pystysiirtymät
    plt.plot(t,u2,'--',t,u4,'--',t,u6,'--')
    plt.subplot(122) # Vaakasiirtymät
    plt.plot(t,u1,'--',t,u3,'--',t,u5,'--')
    plt.show()
    return r

def onedmeshtest6():
    ''' Yksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = F
        
        x0|------x1------|x2
        
    '''    
    
    # Alkutila
    x0 = np.array([0.0, 0.6, 1.0])
    y0 = np.array([0.0, 0.3, 0.0])
    dof = len(x0)*2
    u0 = np.zeros(dof)
    u0[0::2] = x0
    u0[1::2] = y0
    edge_db = [[0,1],[1,2]]

    # Täytyy olla jokin käsitys alkutilasta. Tässä oletetaan että
    # 1) viivaelementit ovat yhtä pitkiä (keskipituus)
    # 2) venymä on poikkeama tästä keskipituudesta

    nedges = len(edge_db)
    
    # Reunan keskipituus
    Ltotal = 0
    i = 0
    for ni,nj in edge_db:
        u0[i] = np.sqrt((x0[ni]-x0[nj])**2 + (y0[ni]-y0[nj])**2)
        Ltotal += u0[i]
        i += 1
    Lavg = Ltotal/nedges
    
    # Kokemus osoittaa, että stabiilissa systeemissä olisi hyvä olla vähän vetoa tai 
    # puristusta, joten säädetään sitä tällä
    k = 0.5
    Lavg *= k

    print("Lavg: {0}".format(Lavg))
    print("x0: {x0}".format(x0=x0))
    print("y0: {y0}".format(y0=y0))
    print("u0: {u0}".format(u0=u0))
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
        v = y[n:2*n]
        dprint("y")
        dprint(y)
        dprint("u")
        dprint(u)
        dprint("v")
        dprint(v)
        dprint("ux")
        dprint(ux)
        dprint("uy")
        dprint(uy)

        # Jäykkyysmatriisi
        def k(edge = None):

            node1 = edge[0]
            node2 = edge[1]

            dprint("node1: {0}    node2: {1}".format(node1,node2))

            x1 = float(x0[node1]) + float(ux[node1])
            y1 = float(y0[node1]) + float(uy[node1])
            x2 = float(x0[node2]) + float(ux[node2])
            y2 = float(y0[node2]) + float(ux[node2])
            
            ux1 = float(ux[node1])
            uy1 = float(uy[node1])
            ux2 = float(ux[node2])
            uy2 = float(ux[node2])
            dprint("uy: {0}".format(uy))
            
            dprint("ux1: {ux1:5.2f}  uy1: {uy1:5.2f}  ux2: {ux2:5.2f}  uy2: {uy2:5.2f}".format(ux1=ux1,uy1=uy1,ux2=ux2,uy2=uy2))
            
            alp = np.arctan2(y2-y1,x2-x1)
            E = 1e6
            A = 1e-6
            EA = 1e3
            L = Lavg#np.sqrt((x2-x1)**2 + (y1-y2)**2)
            
            # "Stress stiffening" kerroin on verrannollinen venymään
            du = np.sqrt((ux2-ux1)**2 + (uy2-uy1)**2)
            dprint("du: {0} , alp: {1} , x1: {x1:5.2f} y1: {y1:5.2f} x2: {x2:5.2f} y2: {y2:5.2f}".format(du,alp,x1=x1,y1=y1,x2=x2,y2=y2))
            
            c = np.cos(alp)
            s = np.sin(alp)
            cc = c*c
            ss = s*s
            cs = c*s
            ki = EA/L*np.matrix([[cc,cs,-cc,-cs],[cs,ss,-cs,-ss],[-cc,-cs,cc,cs],[-cs,-ss,cs,ss]])

#            du = 0
            dprint("ki:")
            dprint(ki*(1+du))
            return ki*(1+du)
        
        def m(edge = None, rho = 1):
            ''' Massamatriisin rakennus '''
            x1 = float(x0[edge[0]])
            y1 = float(y0[edge[0]])
            x2 = float(x0[edge[1]])
            y2 = float(y0[edge[1]])
            L = np.sqrt((x2-x1)**2 + (y1-y2)**2)
            M = rho*L/6*np.matrix([[2,1],[1,2]])
            alp = np.arctan2(y2-y1,x2-x1)
            c = np.cos(alp)
            s = np.sin(alp)
            B = np.matrix([[c,s,0,0],[0,0,c,s]])
            return B.T*M*B

        k1 = k(edge = edge_db[0])
        k2 = k(edge = edge_db[1])
        m1 = m(edge = edge_db[0])
        m2 = m(edge = edge_db[1])
        
        V1 = np.matrix([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0]])
        V2 = np.matrix([[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])

        dprint(V1)
        dprint(k1)
        dprint(V2)
        dprint(k2)

        dof = 6 # Vapausasteiden määrä

        m = np.ones(dof)*1
        dprint(np.matrix(np.diag(m)))

        K = np.matrix(np.zeros((dof,dof)))#V1.T*k1*V1 + V2.T*k2*V2
        dprint("K")
        dprint(K)
        M = np.matrix(np.diag(m)) # + V1.T*m1*V1 + V2.T*m2*V2 
        dprint("M")
        dprint(M)
        
        invM = np.linalg.inv(M)
        # Vaimennusmatriisi
        alpha = 1
        beta = 0
        
        C = alpha*M + beta*K
        
        A = np.matrix(np.zeros((2*dof,2*dof)))
        A[0:dof,dof:2*dof] = np.diag(np.ones(n))
        A[dof:2*dof,0:dof] = -invM*K
        A[dof:2*dof,dof:2*dof] = -invM*C
        dprint("A")
        dprint(A)
        
        def getForceEq(nodeid):
            ''' Palauttaa voimaresultantin solmupisteessä i '''
            x = x0[nodeid] + ux[nodeid]
            y = y0[nodeid] + uy[nodeid]
            print("node: {node}   x: {nodex}   y: {nodey}".format(node = nodeid, nodex = x, nodey = y))
            
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
                dL = L - Lavg
#                print("dL(edge={edge}) = {length}".format(edge=edge_idx,length=dL))
                alp = np.arctan2(yn-y,xn-x)
                print("alpha(edge={edge}) = {alpha}".format(edge=edge_idx,alpha=alp*180/np.pi))
                k = 1
                F = k*dL
                Fx = F*np.cos(alp)
                Fy = F*np.sin(alp)
                print("F(edge={edge}) : Fx,Fy = {Fx:5.2f} , {Fy:5.2f} . Resultant: {Force}".format(edge=edge_idx,Force=F,Fx=Fx,Fy=Fy))
                Fx_tot += Fx
                Fy_tot += Fy
            
            dprint("Fx: {Fx:5.2f}  Fy: {Fy:5.2f}".format(Fx=Fx_tot,Fy=Fy_tot))
            return Fx_tot,Fy_tot
        
        # Voimavektori
        F = np.zeros(dof)
        Fx,Fy = getForceEq(1)
        F[2] = Fx
        F[3] = Fy
        
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
        # Eliminoidaan vapausasteita vaikkapa tällä tavalla
        lock = np.ones(dof*2)
        lock[0] = 0
        lock[1] = 0
        lock[4] = 0
        lock[5] = 0
        dprint("g: {0}".format(g))
        dprint("lock: {0}".format(lock))
        g *= lock
        dprint("g (new): {0}".format(g))
        return g.T
    
#    dt = 0.1
#    t0 = 0
    #y0[0:dof] = u0
##
#    g0 = g(t0,y0)
#    print("g0")
#    print g0
#    g1 = g(t0+dt,g0)
#    print("g1")
#    print g1
   
    r = ode(g)
#    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(np.zeros(2*dof))
##    r.set_f_params([m,c,k])
    t1 = 10
    dt = 0.1
    t = []
    y = []
    condition = True
    i = 0
    itermax = 1000
##    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y.real)
        print r.t, r.y.real
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False


    y = np.array(y)

    x2 = y[:,2] + x0[1]
    y2 = y[:,3] + y0[1]


    import matplotlib.pylab as plt
    plt.subplot(121) # X
    plt.plot(t,x2,'--')
    plt.subplot(122) # Y
    plt.plot(t,y2,'--')
    plt.show()

#    import matplotlib.pylab as plt
#    plt.subplot(121) # X
#    plt.plot(t,x1,'--',t,x2,'--',t,x3,'--')
#    plt.subplot(122) # Y
#    plt.plot(t,y1,'--',t,y2,'--',t,y3,'--')
#    plt.show()

    return r

def onedmeshtest7():
    
    ''' Kaksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = f
        
        kuva prujussa.
        
    '''    
    
    # Alkutila
    x0 = np.array([np.random.rand(), -1.0, 1.0, 0.0])
    y0 = np.array([np.random.rand(), 1.0, 1.0, -1.0])
    dof = len(x0)*2
    u0 = np.zeros(dof)
    u0[0::2] = x0
    u0[1::2] = y0
    edge_db = [[0,1],[0,2],[0,3]]

    # Täytyy olla jokin käsitys alkutilasta. Tässä oletetaan että
    # 1) viivaelementit ovat yhtä pitkiä (keskipituus)
    # 2) venymä on poikkeama tästä keskipituudesta

    nedges = len(edge_db)
    
    # Reunan keskipituus
    Ltotal = 0
    i = 0
    for ni,nj in edge_db:
        u0[i] = np.sqrt((x0[ni]-x0[nj])**2 + (y0[ni]-y0[nj])**2)
        Ltotal += u0[i]
        i += 1
    Lavg = Ltotal/nedges
    
    # Kokemus osoittaa, että stabiilissa systeemissä olisi hyvä olla vähän vetoa tai 
    # puristusta, joten säädetään sitä tällä
    k = 0.5
    Lavg *= k

    print("Lavg: {0}".format(Lavg))
    print("x0: {x0}".format(x0=x0))
    print("y0: {y0}".format(y0=y0))
    print("u0: {u0}".format(u0=u0))
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
        v = y[n:2*n]
        dprint("y")
        dprint(y)
        dprint("u")
        dprint(u)
        dprint("v")
        dprint(v)
        dprint("ux")
        dprint(ux)
        dprint("uy")
        dprint(uy)

        dof = 8 # Vapausasteiden määrä

        m = np.ones(dof)*1
        dprint(np.matrix(np.diag(m)))

        K = np.matrix(np.zeros((dof,dof)))
        dprint("K")
        dprint(K)
        M = np.matrix(np.diag(m)) # + V1.T*m1*V1 + V2.T*m2*V2 
        dprint("M")
        dprint(M)
        
        invM = np.linalg.inv(M)

        # Vaimennusmatriisi
        alpha = 1
        beta = 0
        
        C = alpha*M + beta*K
        
        A = np.matrix(np.zeros((2*dof,2*dof)))
        A[0:dof,dof:2*dof] = np.diag(np.ones(n))
        A[dof:2*dof,0:dof] = -invM*K
        A[dof:2*dof,dof:2*dof] = -invM*C
        dprint("A")
        dprint(A)
        
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
                dL = L - Lavg
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
            
            dprint("Fx: {Fx:5.2f}  Fy: {Fy:5.2f}".format(Fx=Fx_tot,Fy=Fy_tot))
            return Fx_tot,Fy_tot
        
        # Voimavektori
        F = np.zeros(dof)
        Fx,Fy = getForceEq(0)
        F[0] = Fx
        F[1] = Fy
        
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
        # Eliminoidaan vapausasteita vaikkapa tällä tavalla
        lock = np.ones(dof*2)
        lock[2] = 0
        lock[3] = 0
        lock[4] = 0
        lock[5] = 0
        lock[6] = 0
        lock[7] = 0
        dprint("g: {0}".format(g))
        dprint("lock: {0}".format(lock))
        g *= lock
        dprint("g (new): {0}".format(g))
        return g.T
    
   
    r = ode(g)
#    r.set_integrator('zvode', method='bdf')
    r.set_initial_value(np.zeros(2*dof))
##    r.set_f_params([m,c,k])
    t1 = 10
    dt = 0.1
    t = []
    y = []
    condition = True
    i = 0
    itermax = 1000
##    from scipy.integrate import ode
    while r.successful() and condition:
        r.integrate(r.t+dt)
        t.append(r.t)
        y.append(r.y.real)
        print r.t, r.y.real
        i += 1
        if r.t >= t1 or i >= itermax: 
            condition = False


    y = np.array(y)

    xx = y[:,0] + x0[0]
    yy = y[:,1] + y0[0]


    import matplotlib.pylab as plt
    plt.subplot(121) # X
    plt.plot(t,xx,'--')
    plt.subplot(122) # Y
    plt.plot(t,yy,'--')
    plt.show()

#    import matplotlib.pylab as plt
#    plt.subplot(121) # X
#    plt.plot(t,x1,'--',t,x2,'--',t,x3,'--')
#    plt.subplot(122) # Y
#    plt.plot(t,y1,'--',t,y2,'--',t,y3,'--')
#    plt.show()

    return r
