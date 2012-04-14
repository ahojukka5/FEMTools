# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:47:02 2012
@author: Jukka Aho

Pintarakenteet, harjoitus 6, tehtävä 1

"""

from __future__ import division
from teefem import *
from teefem.models.dkt import *
import matplotlib.pylab as plt

t = lambda k,e: 0.6
p0 = 10.0e3
E = 30.0e9
nu = 0.2

mesh0 = mesh(filename = 'geom.msh')
#print mesh.status
mesh0.create_node_groups()
#print mesh.status
mdl = dkt(mesh = mesh0)

OM = mdl.elset('OM1','OM2','OM3','OM4','OM5')

mat = materials.elastic(E = E, nu = nu)
assign_material(elements = OM, material = mat)

# Jakautunut kuormitus
bc1 = pressure_bc(pressure = lambda k,e: -p0)
assign_bc(elements = OM, bc = bc1)

# Jäykkä tuenta pilarit + vasen reuna
bc2 = dirichlet_bc(encastre = True)
bc3 = dirichlet_bc(dz = 0)
assign_bc(nodes = mdl.nset('OM1'), bc = bc2)
assign_bc(nodes = mdl.nset('OM2'), bc = bc2)
assign_bc(nodes = mdl.nset('OM3'), bc = bc2)
assign_bc(nodes = mdl.nset('OM4'), bc = bc2)
assign_bc(nodes = mdl.nset('GA8'), bc = bc2)

# Ylä- ja alareunalla niveltuettu
assign_bc(nodes = mdl.nset('GA5'), bc = bc3)
assign_bc(nodes = mdl.nset('GA7'), bc = bc3)

carel = platechar(thickness = t)
assign_char(elements = OM, char = carel)

mdl.static_solve()

# Jälkikäsittely
#mdl.plot()
#plt.show()

###############################################################################
############################# Jälkikäsittely ##################################
###############################################################################

def plotmdl(mdl):
    ''' Plottailee kaikenlaisia kuvaajia '''

    # Mises jännitys integroimispisteitä interpoloimalla
    # http://www.scipy.org/Cookbook/Matplotlib/Gridding_irregularly_spaced_data
    from scipy.interpolate import griddata
    import numpy as np


    # Primäärikenttä solmupisteissä
    xx = []
    yy = []
    dz = []
    drx = []
    dry = []
    for n in mdl.nodes:
        xx.append(n.x)
        yy.append(n.y)
        dz.append(n.fields['DEPL']['DZ'])
        drx.append(n.fields['DEPL']['DRX'])
        dry.append(n.fields['DEPL']['DRY'])

    # Sekundäärikenttä integroimispisteissä
    x = []
    y = []
    Mx = []
    My = []
    Mxy = []
    Sx = []
    Sy = []
    Sxy = []
    Vmis = []
    for e in OM:
        for ke in e.ipoints:
            x.append(e.geom.x(*ke))
            y.append(e.geom.y(*ke))
            M = e.M(*ke)
            Mx.append(M[0,0])
            My.append(M[1,0])
            Mxy.append(M[2,0])
#            S = e.S(*ke)/1e6
#            Sx.append(S[0,0])
#            Sy.append(S[1,0])
#            Sxy.append(S[1,0])
#            Vmis.append(e.vmis(*ke)/1e6)
    
    # Plottaa kentän ja tallentaa tiedostoon
    def plot(x,y,field,filename,c=200):
        plt.figure()
        # define grid.
        xi = np.linspace(min(x),max(x),100)
        yi = np.linspace(min(y),max(y),100)
        # grid the data.
        si_lin = griddata((x, y), field, (xi[None,:], yi[:,None]), method='linear')
        si_cub = griddata((x, y), field, (xi[None,:], yi[:,None]), method='linear')
        print np.min(field)
        print np.max(field)
        plt.subplot(211)
        # contour the gridded data, plotting dots at the randomly spaced data points.
        CS = plt.contour(xi,yi,si_lin,c,linewidths=0.5,colors='k')
        CS = plt.contourf(xi,yi,si_lin,c,cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar
        # plot data points.
        #    plt.scatter(x,y,marker='o',c='b',s=5)
        plt.xlim(min(x),max(x))
        plt.ylim(min(y),max(y))
        plt.title('Lineaarinen interpolointi')
        #plt.tight_layout()
        plt.subplot(212)
        # contour the gridded data, plotting dots at the randomly spaced data points.
        CS = plt.contour(xi,yi,si_cub,c,linewidths=0.5,colors='k')
        CS = plt.contourf(xi,yi,si_cub,c,cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar
        # plot data points.
        #    plt.scatter(x,y,marker='o',c='b',s=5)
        plt.xlim(min(x),max(x))
        plt.ylim(min(y),max(y))
        plt.title('Kuubinen interpolointi')
        plt.savefig(filename)
        #plt.tight_layout()

#    plot(xx,yy,dz,'output/DZ.pdf')
#    plot(xx,yy,drx,'output/DRX.pdf')
#    plot(xx,yy,dry,'output/DRY.pdf')
    plot(x,y,Mx,'output/Mx.pdf')
    plot(x,y,My,'output/My.pdf')
    plot(x,y,Mxy,'output/Mxy.pdf')
#    plot(x,y,Sx,'output/Sx.pdf')
#    plot(x,y,Sy,'output/Sy.pdf')
#    plot(x,y,Sxy,'output/Sxy.pdf')
#    plot(x,y,Vmis,'output/Vmis.pdf')


plotmdl(mdl)


#
##print mdl.group_no
##print dir(mdl.group_no)
#P1 = iter(mdl.group_no['P1']).next()
#print("DY: %0.14E"%(P1.fields['DEPL']['DY']))