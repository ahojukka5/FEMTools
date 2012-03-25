# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:47:02 2012
@author: Jukka Aho

Tällä voi tutkia eri laattamallien käyttöä pyörähdyssymmetriseen rakenteeseen.

"""

import teefem
import teefem.models


# Lähtötiedot
t = 25e-3
E = 100e9
nu = 0.3
q = 100e3
a = 1
D = E*t**3/(12*(1-nu**2))
model = teefem.models.dkt
#model = teefem.models.mindlin
#model = teefem.models.mindlin_reduced
#model = teefem.models.minth # Tessler & Hughes
#model = teefem.models.dsts6

w = lambda r: -a**4*q*(1 - r**2/a**2)**2/(64*D)
vmis = lambda r: (36*D**2*(a**2*q*(1 - r**2/a**2)/(8*D) - q*r**2/(8*D))**2/t**4 - 36*D**2*(a**2*q*(1 - r**2/a**2)/(8*D) - q*r**2/(8*D))*(a**2*nu*q*(1 - r**2/a**2)/(16*D) + a**2*q*(1 - r**2/a**2)/(16*D) - q*r**2/(8*D))/t**4 + 36*D**2*(a**2*nu*q*(1 - r**2/a**2)/(16*D) + a**2*q*(1 - r**2/a**2)/(16*D) - q*r**2/(8*D))**2/t**4)**(1/2)
mesh = teefem.mesh.unitcircle(R=1)
mat = teefem.materials.elastic(E = E, nu = nu)
mdl = model(mesh = mesh)

OM_el = mdl.elset['OM1']
GA1_no = mdl.nset['GA1']

teefem.assign_material(elements = OM_el, material = mat)

bc1 = teefem.pressure_bc(pressure = lambda k,e: -q)
bc2 = teefem.dirichlet_bc(encastre = True)

teefem.assign_bc(elements = OM_el, bc = bc1)
teefem.assign_bc(nodes = GA1_no, bc = bc2)

carel = teefem.plate_functions.PlateCharacteristic(thickness = lambda k,e: t)
teefem.assign_char(elements = OM_el, char = carel)

mdl.static_solve()

dymin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])

print("DZ (acc)  : %0.14E"%(w(0)))
print("DZ        : %0.14E"%(dymin))

import matplotlib.pylab as plt
mdl.plot()
#plotmdl(mdl, vmis)
plt.show()