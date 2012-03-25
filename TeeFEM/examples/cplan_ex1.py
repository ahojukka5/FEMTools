# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:47:02 2012
@author: Jukka Aho

Levyongelma

         q(x) = q0
   ||||||||||||||||||||||||
   --------- GA2 ----------
   |                      |
   |                      |
  GA1        OM1          |
   |                      |
   |                      |
   ----------------------P1
                          |
                          |F

"""

import teefem
import numpy as np

#meshfile = teefem.os.path.join(teefem.datadir, 'cplan_tria3.mail') # Linear triangle elements
meshfile = teefem.os.path.join(teefem.datadir, 'cplan_tria6.mail') # Quadratic triangle elements
#meshfile = teefem.os.path.join(teefem.datadir, 'cplan_quad4.mail') # Linear quadrangle elements
#meshfile = teefem.os.path.join(teefem.datadir, 'cplan_quad8.mail') # Quadratic quadrangle elements

t = 10e-3
F = 100e3
q = 100e3

mesh = teefem.geom.mesh(filename = meshfile)

mdl = teefem.models.cplan(mesh = mesh)

# Ryhmät
OM1_el = mdl.elset['OM1'] # Koko domain
GA1_no = mdl.nset['GA1'] # Vasen reuna
GA2_el = mdl.elset['GA2'] # Yläreuna
P1_no = mdl.nset['P1'] # Pistevoima

# Materiaali
mat = teefem.materials.elastic(E = 210.0e9, nu = 0.3)
teefem.assign_material(elements = OM1_el, material = mat)

# Reunaehdot, huom. TeeFEMissä tasojännitystila-elementin paksuus on oletuksena 1
bc1 = teefem.pressure_bc(pressure = lambda k: -q/t)
bc2 = teefem.dirichlet_bc(encastre = True)
bc3 = teefem.nodal_force(fy = -F/t)

# Kytketään reunaehdot malliin
# 1. Jakautunut kuormitus elementtiryhmään GA2
teefem.assign_bc(elements = GA2_el, bc = bc1)
# 2. Vasen reuna jäykästi kiinni
teefem.assign_bc(nodes = GA1_no, bc = bc2)
# 3. Pistevoima alareunaan
teefem.assign_bc(nodes = P1_no, bc = bc3)

# Ratkaisu
import time
t0 = time.clock()
mdl.static_solve()
t1 = time.clock()
print("Solving KU=F took {0} seconds".format(t1-t0))

######################
### Jälkikäsittely ###
######################

# Taipuma solmupisteryhmässä P1
for node in P1_no:
    print("DY: %0.14E"%(node.fields['DEPL']['DY']))

sf = 50

import matplotlib.pylab as plt

# Siirtymäkenttä solmupisteissä matplotlibillä

plt.figure(1)
x = [n.x+n.fields['DEPL']['DX']*sf for n in mdl.nodes]
y = [n.y+n.fields['DEPL']['DY']*sf for n in mdl.nodes]
plt.scatter(x,y)

# Mises jännitys integroimispisteitä interpoloimalla
# http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
plt.figure(2)
from scipy.interpolate import griddata
x = []
y = []
vmis = []
for e in mdl.elset['OM1']:
    for ke in e.ipoints:
        x.append(e.geom.x(*ke))
        y.append(e.geom.y(*ke))
        vmis.append(e.vmis(*ke)/1e6)
# define grid.
xi = np.linspace(min(x),max(x),100)
yi = np.linspace(min(y),max(y),100)
# grid the data.
si = griddata((x, y), vmis, (xi[None,:], yi[:,None]), method='cubic')
# contour the gridded data, plotting dots at the randomly spaced data points.
CS = plt.contour(xi,yi,si,30,linewidths=0.5,colors='k')
CS = plt.contourf(xi,yi,si,30,cmap=plt.cm.jet)
plt.colorbar() # draw colorbar
# plot data points.
# plt.scatter(x,y,marker='o',c='b',s=5)
plt.xlim(min(x),max(x))
plt.ylim(min(y),max(y))
#plt.tight_layout()
plt.show()
