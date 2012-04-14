# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 19:10:42 2012

@author: Jukka Aho

Rotsysolid testi

Yla ja alareuna, GA1
Keskipiste oikealla P1
Domain OM1

"""

from teefem import *
from teefem.models.rotsysolid import *
from numpy import *
from matplotlib.pylab import *
import teefem

# Lähtötiedot
t = 0.010
E = 210.0e9
nu = 0.3
F = 10e3

# Ladataan verkko
mesh = mesh(filename = 'rotsysolid_test.msh')

# Solmupisteryhmät
mesh.create_node_groups()

# Malli
mdl = ROTSOL(mesh = mesh)

# Materiaali
assign_material(
    elements = mdl.elset('OM1'),
    material = materials.elastic(E = E, nu = nu),
    )

# Pistekuorma keskelle
assign_bc(
    nodes = mdl.nset('P1'), 
    bc = teefem.nodal_force(fx = -F),
    )

# Jäykkä tuki ylä- ja alareunaan
assign_bc(
    nodes = mdl.nset('GA1'),
    bc = dirichlet_bc(dx=0,dy=0),
    )

# Ratkaistaan maksimitaipumaa alustakertoimen funktiona.
# Saadaan kiva animaatio. Ja kuvaaja.

mdl.static_solve()
mdl.plot(scalefactor=1e6)
show()