# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 04:56:18 2012

@author: Jukka Aho

Geometrinen epälineaarisuus. Jousi, joka on toisesta päästä niveltuettu. Esim.
vetokoe. 

Ratkaisee staattisen tasapainon k(u)u = f (stress stiffening??)

Esimerkiksi k(u) = E*A/L*u ( ei mitään fysikaalista merkitystä )

Tästä voisi osaava tehdä version jossa olisi materiaalinen epälineaarisuus tmv.

"""

import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import newton

E = 210e9 # Materiaali
d = 15e-3 # Sauvan poikkipinta-alan halkaisija
A = np.pi*d**2/4 # Sauvan poikkipinta-ala
L = 1 # Sauvan pituus

FF = 100e3 # Maksimivoima
n = 30 # Kuormitusaskelet

def sol(k):
    F = np.linspace(0,FF,n)
    F_ = []
    u_ = []
    eps_ = []
    sig_ = []
    u = 0
    i = 1
    for Fi in F:
        func = lambda u: k(u)*u-Fi
        usol = newton(func,u)
        print("Sol {0:2d}: F = {1:0.2f}, u = {2:0.2f}".format(i,Fi,usol))
        F_.append(Fi)
        u_.append(usol)
        eps_.append(usol/L*1e6)
        sig_.append(Fi/A/1e6)
        i += 1
    return F_,u_,eps_,sig_

def plot(F,u,eps,sig):
    plt.figure()
    plt.subplot(121)
    plt.plot(F,u,'--o')
    plt.xlabel('Force')
    plt.ylabel('Displacement')
    plt.subplot(122)
    plt.plot(eps,sig,'--o')
    plt.xlabel('Strain [u]')
    plt.ylabel('Stress [MPa]')

F1,u1,eps1,sig1 = sol(lambda u: E*A/L) # Jäykkyys venymän funktiona, lineaariteoria

F2,u2,eps2,sig2 = sol(lambda u: u*E*A/L) # Jäykkyys venymän funktiona, lineaariteoria

plot(F1,u1,eps1,sig1)
plot(F2,u2,eps2,sig2)

plt.show()