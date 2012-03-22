# -*- coding: utf-8 -*-
"""
Created on Fri Mar 16 13:57:08 2012

@author: Jukka

Reuna-arvotehtävä:
    u'' + u = f     f e R
    u(x0) = u0
    u'(x1) = u1

Tarkka ratkaisu:
    u(x) = c1*cos(x) + c2*sin(x) + f
    Ac = b
    A = [[cos(x0), sin(x0)],[-sin(x1),cos(x1)]]
    c = [c1,c2].T
    b = [u0, u1].T

Ratkaisee elementtimenetelmällä. Doplhin ratkaisu ei toimi.

"""


from __future__ import division

f = 1.0
x0 = 0.0
x1 = 1.0
u0 = 0.0
u1 = 1.0

def femsolve():
    
    ''' Bilineaarinen muoto:

        a(u,v) = L(v)
        a(u,v) = (inner(grad(u), grad(v)) + u*v)*dx
        L(v) = f*v*dx - g*v*ds
        g(x) = -du/dx = -u1, x = x1
        u(x0) = u0
        Omega = {xeR|x0<=x<=x1}

    '''

    from dolfin import UnitInterval, FunctionSpace, DirichletBC, TrialFunction
    from dolfin import TestFunction, grad, Constant, Function, solve, inner, dx, ds
    from dolfin import MeshFunction, assemble
    import dolfin
#    from dolfin import set_log_level, PROCESS

    # Create mesh and define function space
    mesh = UnitInterval(30)
    V = FunctionSpace(mesh, 'Lagrange', 2)

    boundaries  = MeshFunction('uint', mesh, mesh.topology().dim()-1)

    boundaries.set_all(0)

    class Left(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            tol = 1E-14   # tolerance for coordinate comparisons
            return on_boundary and abs(x[0]) < tol

    class Right(dolfin.SubDomain):
        def inside(self, x, on_boundary):
            return dolfin.near(x[0], 1.0)
    
    left = Left()
    right = Right()
    
    left.mark(boundaries, 1)
    right.mark(boundaries, 2)

#    def u0_boundary(x):
#        return abs(x[0]) < tol
#    
#    bc = DirichletBC(V, Constant(u0), lambda x: abs(x[0]) < tol)
    
    bcs = [DirichletBC(V, Constant(u0), boundaries, 1)]
    
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    a = (inner(grad(u), grad(v)) + u*v)*dx
    g = Constant(-u1)
    L = Constant(f)*v*dx - g*v*ds(2)
    
 #   set_log_level(PROCESS)
    # Compute solution
    
    A = assemble(a, exterior_facet_domains=boundaries)
    b = assemble(L, exterior_facet_domains=boundaries)
    for bc in bcs: 
        bc.apply(A, b)
    
    u = Function(V)
    solve(A, u.vector(), b, 'lu')
    
    coor = mesh.coordinates()
    u_array = u.vector().array()
    a = []
    b = []
    for i in range(mesh.num_vertices()):
        a.append(coor[i])
        b.append(u_array[i])
        print('u(%3.2f) = %0.14E'%(coor[i],u_array[i]))
    
    import numpy as np
    np.savez('fem',a,b)
#    print np.array(coor)
#    print np.array(u_array)

    # Verification
#    u_exact = Expression('x[0]*x[0]')
#    u_e = interpolate(u_exact, V)
#    print 'Max error:', \
#          numpy.abs(u_e.vector().array() - u.vector().array()).max()
#    
#    
#    plot(u)
#    plot(mesh)
#    interactive()

def accsolve():
    import numpy as np
    A = np.array([[np.cos(x0), np.sin(x0)],[-np.sin(x1),np.cos(x1)]])
    b = np.array([u0-f, u1])
    c = np.linalg.solve(A,b)
    print c
    u_acc = np.vectorize(lambda x: c[0]*np.cos(x) + c[1]*np.sin(x) + f)
    du_acc = np.vectorize(lambda x: -c[0]*np.sin(x) + c[1]*np.cos(x))
    print("du_acc({0}) = {1}".format(x1, du_acc(x1)))
    x = np.linspace(x0,x1)
    u = u_acc(x)
    np.savez('acc',x,u)

def plot():
    import numpy as np
    acc = np.load('acc.npz')
    fem = np.load('fem2.npz')
    acc_x = acc['arr_0']
    acc_u = acc['arr_1']
    fem_x = fem['arr_0']
    fem_u = fem['arr_1']
#    print acc_x
#    print acc_u
    print fem_x
    print fem_u
    import matplotlib.pylab as plt
    plt.plot(acc_x, acc_u,'b--', fem_x,fem_u,'ro')
    plt.show()

def fem2():
    import numpy as np
    n = 5
    h = (x1-x0)/n
    k = lambda h: h/6*np.array([[2,1],[1,2]]) - 1/h*np.array([[1,-1,],[-1,1]]) 
    fq = lambda h: f*h/2*np.array([1,1])
    dim = n+2
    K = np.zeros((dim,dim))
    R = np.zeros(dim)
    for i in range(n):
        start = i
        end = i+2
        K[start:end,start:end] += k(h)
        R[start:end] += fq(h)
    # Reunaehdot
    K[-1,0] = 1
    K[0,-1] = 1

    # Kuormitusvektori
    R[-1] += u0
    R[-2] += -u1

    sol = np.linalg.solve(K,R)
    y = sol[0:n+1]
    x = np.linspace(x0,x1,len(y))
    np.savez('fem2',x,y)
#    import matplotlib.pylab as plt
#    plt.plot(x,y)
#    plt.show()

accsolve()
#femsolve()
fem2()
plot()