# -*- coding: utf-8 -*-
"""
Created on Fri Mar 23 18:43:43 2012

@author: Jukka Aho

Esimerkki "stress stiffening" josta voi rakentaa dynamiikka solverin.

"""

import numpy as np
from numpy import matrix,array

def onedmeshtest5():
    ''' Yksidimensioinen verkon tasapainotustesti 
        mu'' + cu' + ku' = F
        
        x0|------x1------|x2
        
    '''    
    
    # Alkutila
    L = 1
    x0 = np.array([0,0.5,L])
    y0 = np.array([0,0,0])
    edge_db = [[0,1],[1,2]]
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
            E = 210e9
            A = 1e-6
            L = 0.5#np.sqrt((x2-x1)**2 + (y1-y2)**2)
            
            # "Stress stiffening" kerroin on verrannollinen venymään
            du = np.sqrt((ux2-ux1)**2 + (uy2-uy1)**2)
            dprint("du: {0} , alp: {1} , x1: {x1:5.2f} y1: {y1:5.2f} x2: {x2:5.2f} y2: {y2:5.2f}".format(du,alp,x1=x1,y1=y1,x2=x2,y2=y2))
            
            c = np.cos(alp)
            s = np.sin(alp)
            cc = c*c
            ss = s*s
            cs = c*s
            ki = E*A/L*np.matrix([[cc,cs,-cc,-cs],[cs,ss,-cs,-ss],[-cc,-cs,cc,cs],[-cs,-ss,cs,ss]])

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
        a1 = matrix([1,0,0,0,0,0])
        a2 = matrix([0,1,0,0,0,0])
        a3 = matrix([0,0,0,0,1,0])
        a4 = matrix([0,0,0,0,0,1])
        s = 1e6
        K += s*(a1.T*a1 + a2.T*a2 + a3.T*a3 + a4.T*a4)

        
        
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
        g = A*yy + b
        #print g
        return g
    
#    dt = 0.1
#    t0 = 0
    y0 = np.zeros(12)
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
    plt.plot(t,u2,'--',t,u4,'--',t,u6,'--')
    plt.show()
    return r

if __name__ == '__main__':
    #unitsquare()
    #sgmtest()
    #trimeshtest()
    onedmeshtest5()
