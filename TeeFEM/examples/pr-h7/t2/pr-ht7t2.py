# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 23:15:29 2012

@author: Jukka Aho
"""

from __future__ import division
import teefem
from teefem.models import dkt, cplan
import numpy as np
import matplotlib.pylab as plt

t = 10e-3
E = 210e9
nu = 0.3
a = 1
b = 1
q = 1e6

def malli0():
    ''' Tasojännitystilan ratkaisu. '''
    mesh = teefem.geom.mesh(filename = 'cplangeo.msh')
    mdl = cplan.cplan(mesh = mesh)
    OM1_el = mdl.elset['OM1']
    mat = teefem.materials.elastic(E = 210.0e9, nu = 0.3)
    teefem.assign_material(elements = OM1_el, material = mat)
#    nset1 = set().union(mdl.nset['GA1'], mdl.nset['GA2'], mdl.nset['GA3'], mdl.nset['GA4'])
#    nset1 = set().union(mdl.nset['GA2'], mdl.nset['GA4'])
    teefem.assign_bc(
        nodes = mdl.nset['GA5'], 
        bc = teefem.dirichlet_bc(dx=0,dy=0),
        )
    teefem.assign_bc(
        elements = mdl.elset['GA1'],
        bc = teefem.pressure_bc(pressure = lambda k: +q), # TJT oletuspaksuus t=1
        )
    teefem.assign_bc(
        elements = mdl.elset['GA3'],
        bc = teefem.pressure_bc(pressure = lambda k: -q), # TJT oletuspaksuus t=1
        )
    mdl.static_solve()
    
    ### Jälkikäsittelyä ###
    
    from scipy.interpolate import griddata
    
    # Primäärikenttä solmupisteissä
    x = []
    y = []
    dx = []
    dy = []
    for n in mdl.nodes:
        x.append(n.x)
        y.append(n.y)
        dx.append(n.fields['DEPL']['DX'])
        dy.append(n.fields['DEPL']['DY'])

    xx = []
    yy = []
    Nx = []
    Ny = []
    Nxy = []
    vmis = []
    for e in mdl.elset['OM1']:
        for ke in e.ipoints:
            xx.append(e.geom.x(*ke))
            yy.append(e.geom.y(*ke))
            S = e.S(*ke)
            Nx.append(S[0,0])
            Ny.append(S[1,0])
            Nxy.append(S[2,0])
            vmis.append(e.vmis(*ke))

    def plot(x,y,data,filename):
        plt.figure()
        # define grid.
        xi = np.linspace(min(x),max(x),100)
        yi = np.linspace(min(y),max(y),100)
        # grid the data.
        si = griddata((x, y), data, (xi[None,:], yi[:,None]), method='linear')
        # contour the gridded data, plotting dots at the randomly spaced data points.
        CS = plt.contour(xi,yi,si,30,linewidths=0.5,colors='k')
        CS = plt.contourf(xi,yi,si,30,cmap=plt.cm.jet)
        plt.colorbar() # draw colorbar
        # plot data points.
        # plt.scatter(x,y,marker='o',c='b',s=5)
        plt.xlim(min(x),max(x))
        plt.ylim(min(y),max(y))
        plt.savefig(filename)
        #plt.tight_layout()
    
    plot(x,y,dx,'dx.pdf')
    plot(x,y,dy,'dy.pdf')
    
    plot(xx,yy,Nx,'Nx.pdf')
    plot(xx,yy,Ny,'Ny.pdf')
    plot(xx,yy,Nxy,'Nxy.pdf')
    
    plt.show()

        
def malli1():

    ''' DKT suorakaidelaatan lommahdus. '''


    mesh = teefem.geom.mesh(filename = 'geom.msh')
    mat = teefem.materials.elastic(E = E, nu = nu)
    mdl = dkt.dkt(mesh = mesh)
    
    OM_el = mdl.elset['OM1']

    teefem.assign_material(elements = OM_el, material = mat)
    
    nset1 = set().union(mdl.nset['GA1'], mdl.nset['GA2'], mdl.nset['GA3'], mdl.nset['GA4'])
    teefem.assign_bc(
        nodes = nset1, 
        bc = teefem.dirichlet_bc(dz = 0),
        )

    carel = teefem.plate_functions.platechar(
        thickness = lambda k,e: t,
        Tx = lambda k,e: 1e6,
        )

    teefem.assign_char(elements = OM_el, char = carel)

    mdl.buckling_solve(n_modes = 10)
    for i in range(4):
        mdl.plot_buckling(shapeidx = i)
    
    # Tarkka ratkaisu
    D = E*t**3/(12*(1-nu**2))
    la = lambda m,n: np.pi**2*D/b**2*(m/(a/b) + n**2/m*(a/b))**2
    print("Acc")
    r = range(1,5)
    d = np.array([la(i,j) for i in r for j in r])
    for di in np.sort(d)/1e6:
        print di

    import matplotlib.pylab as plt
    plt.show()

    
#    dymin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])
#    
#    print("DZ (acc)  : %0.14E"%(-w(0)))
#    print("DZ        : %0.14E"%(dymin))
#    print("DZ (CA)   : %0.14E"%(-1.74644966293971E-01))
    
#    import matplotlib.pylab as plt
#    mdl.plot()
#    plotmdl(mdl, vmis)
#    plt.show()

if __name__ == '__main__':
    malli0()