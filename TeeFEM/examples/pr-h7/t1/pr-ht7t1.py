# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 20:31:02 2012

Pintarakenteet, harjoitus 7, tehtävä 1. Moodianalyysi.

@author: Jukka Aho
"""

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from teefem.models.membrane import membrane, memchar
from teefem.geom import mesh
from teefem import assign_bc, dirichlet_bc, assign_char

# Lähtötiedot
m = 1
T = 10
a = 2
b = 1

mesh = mesh(filename = 'geom.msh')

mdl = membrane(mesh = mesh)

# Kaikki reunat kiinni
nset1 = set().union(mdl.nset['GA1'], mdl.nset['GA2'], mdl.nset['GA3'], mdl.nset['GA4'])

assign_bc(
    nodes = nset1, 
    bc = dirichlet_bc(encastre = True)
    )

# Esijännitys
Tf = lambda k,e: T
car = memchar(Tx = Tf, Ty = Tf, m = lambda k,e: m)
assign_char(elements = mdl.elements, char = car)

mdl.modal_solve(n_modes = 4, export_matrices = False)

for i in range(4):
    mdl.plot_modal(shapeidx = i)

# Tarkka ratkaisu
w = lambda i,j: np.sqrt(T/m)*np.sqrt((i*np.pi/a)**2 + (j*np.pi/b)**2)
print("acc:")
X = range(1,4)
for i in X:
    for j in X:
        print("i: {0}  j: {1}  w: {2}".format(i,j,w(i,j)))

plt.show()