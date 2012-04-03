# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 19:10:42 2012

@author: Jukka Aho

Pintarakenteet, harjoitus 8, tehtävä 1 b)

Jäykistetty laatta (ortotrooppinen materiaali)

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
B = 0.200
H = 0.050

# Ladataan verkko
mesh = mesh(filename = 'prht8t1b.msh')

# Solmupisteryhmät
mesh.create_node_groups()

# Malli
mdl = dkt(mesh = mesh)

# Määritetään ortotrooppinen materiaalimalli
#mat1 = materials.elastic_orthotropic(Ex = Ex, Ey = Ey, nux = nux, nuy = nuy)
mat1 = plate_structural_orthotropic(
    basematerial = materials.elastic(E = E, nu = nu),
    b = t, # Jäykisteen leveys
    B = B, # Jäykistejako
    t = t, # Laatan paksuus
    T = t+H, # Laatan paksuus jäykisteen kohdalta
    )

# Materiaali
assign_material(
    elements = mdl.elset('OM1'),
    material = mat1,
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

# Laatan ominaisuudet
carel = platechar(
    thickness = lambda k,e: t,
#    winkler_c = c,
    )

# Kytketään laatan ominaisuudet laattaelementteihin
assign_char(
    elements = mdl.elset('OM1'),
    char = carel)

# Ratkaisu
mdl.static_solve()

# Jälkikäsittely, taipuma
dzmin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])*1000
print("DZ min: {0:5.4f} mm".format(dzmin))
fig, ax = mdl.plot()
fig.savefig('output/ortotrooppi.png')
show()