# -*- coding: utf-8 -*-
"""
Created on Sat Mar 17 14:31:02 2012

@author: Jukka

Reuna-arvotehtävä:
    -u'' = cos(u)
    u(x0) = u0
    u(pi) = u1

Tarkka ratkaisu: ...?

FEM-ratkaisu ja FDM-ratkaisu

"""

from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import time

x0 = 0
x1 = np.pi
u0 = 0
u1 = 1
itermax = 100
n = 50 # Elementtien / vapausasteiden määrä

def solveFDM():
    ''' Differenssimenetelmä, kiintopisteiteraatio '''
    import numpy as np
#    A = np.matrix(np.diag(2*np.ones(n-1)) + np.diag(-1*np.ones(n-2),k=-1) + np.diag(-1*np.ones(n-2),k=1))
    A = np.zeros((n,n))
    for j in range(1,n-1):
        A[j,j-1] = -1
        A[j,j] = 2
        A[j,j+1] = -1
    A[0,0] = 1 # u0 kerroin
    A[n-1,n-1] = 1 # u1 kerroin
    h = (x1-x0)/n
    u = {}
    u[0] = np.linspace(u0,u1,n)
    x = np.linspace(x0,x1,n)
    i = 0
    while i < itermax:
        ii = i
        print("iter: {0}".format(ii))
        b = np.zeros(n)
        b[0] = u0
        b[n-1] = u1
        for jj in range(1,n-1):
#            print("b[{0}]: {1}".format(jj,b[jj]))
#            print("np.cos(u[{0}][{1}]): {2}".format(ii,jj,np.cos(u[ii][jj])))
            b[jj] = h**2*np.cos(u[ii][jj])
        u[ii+1] = np.linalg.solve(A,b)
#        print u[ii+1]
        i += 1
    
    for ui in u.itervalues():
        plt.plot(x,ui,'--')
    plt.plot(x,u[0],'--o')
    plt.plot(x,u[itermax],'--o')
#    plt.plot(x,u[itermax-1],'--o')
#    plt.plot(x,u[itermax-2],'--o')
#    print("Done")

def solveFEM():
    ''' Elementtimenetelmä, Picardin iteraatio, Gaussin integrointi '''
    import numpy as np
    h = (x1-x0)/n
    ki = 1/h*np.array([[1,-1],[-1,1]])
#    ri = 0.55*h/2*np.array([1,1])
    dim = n+3

    K = np.zeros((dim,dim))
    for i in range(n):
        start = i
        end = i+2
        K[start:end,start:end] += ki
    K[-2,0] = K[0,-2] = 1
    K[-1,-3] = K[-3,-1] = 1

    u = {}
    u[0] = np.linspace(u0,u1,n+1)
    x = np.linspace(x0,x1,n+1)
    
    N = lambda k: np.array([(1-k)/2, (1+k)/2])

    def quad(I,n=2):
#        return I(-1/np.sqrt(3)) + I(1/np.sqrt(3))
        return 2.0*I(0.0)
#        W = (1.0,1.0) # Painokertoimet
#        i = (-1.0/np.sqrt(3.0), 1.0/np.sqrt(3.0)) # Integroimispisteet
#        W = [2.0]
#        P = [0.0]
#        return np.sum([W[i]*I(P[i]) i in range(1)])

    from scipy.integrate import romberg

    def ri(u):
        ''' Epälineaarisen osan integrointi Gaussin kaavalla '''
        J = h/2 # Jakobiaani
        I = lambda k: N(k)*np.cos(np.sum(N(k)*u))*J # Integrandi
        I0 = lambda k: I(k)[0]
        I1 = lambda k: I(k)[1]
        #return 2.0*I(0.0)
        #return quad(I,-1,1)[0]
        return np.array([romberg(I0,-1,1),romberg(I1,-1,1)])
    
    i = 0
    while i < itermax:
        R = np.zeros(dim)
        R[-2] = 0
        R[-1] = 1
        for j in range(n):
            start = j
            end = j+2
            R[start:end] += ri(u[i][start:end])
        u[i+1] = np.linalg.solve(K,R)
        i += 1

    for ui in u.itervalues():
        plt.plot(x,ui[0:n+1],'--')
    plt.plot(x,u[0][0:n+1],'--o')
    plt.plot(x,u[itermax][0:n+1],'--o')
#    print u[itermax][0:n+1]

#    plt.plot(x,y,'--')
#    plt.show()
plt.subplot(121)
solveFDM()
plt.subplot(122)
solveFEM()
plt.tight_layout()
plt.show()
