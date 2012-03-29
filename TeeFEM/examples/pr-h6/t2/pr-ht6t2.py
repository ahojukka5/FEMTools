# -*- coding: utf-8 -*-
"""

Pintarakenteet, harjoitus 6, tehtävä 2

"""

from __future__ import division
import teefem as tf
import matplotlib.pylab as plt

# Lähtötiedot
p0 = 10.0e3
E = 30.0e9
nu = 0.2

# Laatan paksuus
t = lambda k,e: 0.2 # Ohut
#t = lambda k,e: 2.25 # Paksu laatta

# Malli
model = tf.models.dkt
#model = tf.models.dsts6
#model = tf.models.min
#model = tf.models.min_r
#model = tf.models.minth


# Verkko
mesh = tf.geom.mesh(filename = 'tria3_50.msh')
#mesh = tf.geom.mesh(filename = 'tria3_5000.msh')

mesh.create_node_groups()

mdl = model(mesh = mesh)

tf.assign_material(
    elements = mdl.elset['OM1'],
    material = tf.materials.elastic(E = E, nu = nu),
    )

# Jakautunut kuormitus
tf.assign_bc(
    elements = mdl.elset['OM1'], 
    bc = tf.pressure_bc(pressure = lambda k,e: -p0),
    )

# Jäykkä tuenta vasemmalle reunalle (GA4)
tf.assign_bc(
    nodes = mdl.nset['GA4'], 
    bc = tf.dirichlet_bc(encastre = True),
    )

# Niveltuenta ylös ja alas

bcnod = set().union(mdl.nset['GA1'],mdl.nset['GA3'])
tf.assign_bc(nodes = bcnod, bc = tf.dirichlet_bc(dz = 0))

carel = tf.plate_functions.PlateCharacteristic(thickness = t)
tf.assign_char(elements = mdl.elset['OM1'], char = carel)

mdl.static_solve()

# Jälkikäsittely

# Maksimitaipuma
dzmin = min([node.fields['DEPL']['DZ'] for node in mdl.nodes])
print("DZ min: {0:5.4e}".format(dzmin))

#mdl.plot()
plt.show()