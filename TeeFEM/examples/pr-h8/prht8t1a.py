# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 19:10:42 2012

@author: Jukka Aho

Pintarakenteet, harjoitus 8, tehtävä 1 a)

Laatta kimmoisalla alustalla

"""

from teefem import *
from teefem.models.dkt import *
from numpy import *
from matplotlib.pylab import *

# Lähtötiedot
t = 0.010
E = 210.0e9
nu = 0.3
c = 10.0e6
q = 100e3

# Ladataan verkko
mesh = mesh(filename = 'prht8t1a.msh')

# Solmupisteryhmät
mesh.create_node_groups()

# Malli
mdl = dkt(mesh = mesh)

# Materiaali
assign_material(
    elements = mdl.elset('OM1'),
    material = materials.elastic(E = E, nu = nu),
    )

# Jakautunut kuormitus
assign_bc(
    elements = mdl.elset('OM1'),
    bc = pressure_bc(pressure = lambda k,e: -q),
    )

# Niveltuenta reunoille
assign_bc(
    nodes = mdl.nset('GA1', 'GA2', 'GA3','GA4'),
    bc = dirichlet_bc(dz = 0),
    )

# Ratkaistaan maksimitaipumaa alustakertoimen funktiona.
# Saadaan kiva animaatio. Ja kuvaaja.

c_ = []
dzmin_ = []
i = 0

for ki in np.linspace(0,10,100):
    winkler_c = ki*c
    carel = platechar(
        thickness = lambda k,e: t,
        winkler_c = winkler_c)
    assign_char(
        elements = mdl.elset('OM1'),
        char = carel)
    mdl.static_solve()
    dzmin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])*1000
    print("DZ min, c = {0:5.2f} : {1:5.4f} mm".format(winkler_c, dzmin))
    c_.append(winkler_c/1e6)
    dzmin_.append(dzmin)
    fig,ax = mdl.plot()
    fig.savefig('output/winkler_%d.png'%(i))
    i += 1
figure()
plot(c_,dzmin_)
title('Maksimitaipuma alustakertoimen funktiona')
xlabel('Alustakerroin [MPa/m]')
ylabel('Maksimitaipuma [mm]')
savefig('output/taipuma.pdf')