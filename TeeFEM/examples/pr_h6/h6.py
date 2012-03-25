# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:47:02 2012
@author: Jukka Aho

Pintarakenteet, harjoitus 6, tehtävä 1

"""

from __future__ import division
import teefem as tf

t = 0.6
p0 = 10.0e3
E = 30.0e9
nu = 0.2

mesh = tf.geom.mesh(filename = 'geom.msh')
#print mesh.status
mesh.create_node_groups()
#print mesh.status
mdl = tf.models.dkt(mesh = mesh)

#OM = mdl.elements
#OM = mdl.elset['OM1'].union(mdl.elset['OM2'], mdl.elset['OM3'], mdl.elset['OM4'], mdl.elset['OM5'])
OM1 = mdl.elset['OM1']
OM2 = mdl.elset['OM2']
OM3 = mdl.elset['OM3']
OM4 = mdl.elset['OM4']
OM5 = mdl.elset['OM5']
OM = set().union(OM1,OM2,OM3,OM4,OM5)
mdl.elset['OM'] = OM

GA1 = mdl.nset['GA1']
GA2 = mdl.nset['GA2']
GA3 = mdl.nset['GA3']
GA4 = mdl.nset['GA4']
GA5 = mdl.nset['GA5']
GA6 = mdl.nset['GA6']
GA7 = mdl.nset['GA7']
GA8 = mdl.nset['GA8']

mat = tf.materials.elastic(E = E, nu = nu)
tf.assign_material(elements = OM, material = mat)

# Jakautunut kuormitus
bc1 = tf.boundary_conditions.PressureBoundaryCondition(pressure = lambda k,e: -p0)
tf.assign_bc(elements = OM, bc = bc1)

# Jäykkä tuenta pilarit + vasen reuna
bc2 = tf.boundary_conditions.DirichletBoundaryCondition(encastre = True)
bc3 = tf.boundary_conditions.DirichletBoundaryCondition(dz = 0)
tf.assign_bc(nodes = mdl.nset['OM1'], bc = bc2)
tf.assign_bc(nodes = mdl.nset['OM2'], bc = bc2)
tf.assign_bc(nodes = mdl.nset['OM3'], bc = bc2)
tf.assign_bc(nodes = mdl.nset['OM4'], bc = bc2)
tf.assign_bc(nodes = mdl.nset['GA8'], bc = bc2)

# Ylä- ja alareunalla niveltuettu
tf.assign_bc(nodes = mdl.nset['GA5'], bc = bc3)
tf.assign_bc(nodes = mdl.nset['GA7'], bc = bc3)

carel = tf.plate_functions.PlateCharacteristic(thickness = lambda k,e: t)
tf.assign_char(elements = OM, char = carel)

mdl.static_solve()

def plotmdl(mdl):
    ''' Plottailee kaikenlaisia kuvaajia '''

    # Mises jännitys integroimispisteitä interpoloimalla
    # http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata
    import numpy as np

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
    for e in mdl.elset['OM']:
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

    # define grid.
    xi = np.linspace(min(x),max(x),100)
    yi = np.linspace(min(y),max(y),100)
    # grid the data.
    si = griddata((x, y), Vmis, (xi[None,:], yi[:,None]), method='linear')
    print np.min(Vmis)
    print np.max(Vmis)
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,si,30,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,si,30,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
#    plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(min(x),max(x))
    plt.ylim(min(y),max(y))
    #plt.tight_layout()
    plt.show()

plotmdl(mdl)

#import matplotlib.pylab as plt
#mdl.plot()
#plt.show()
#
##print mdl.group_no
##print dir(mdl.group_no)
#P1 = iter(mdl.group_no['P1']).next()
#print("DY: %0.14E"%(P1.fields['DEPL']['DY']))