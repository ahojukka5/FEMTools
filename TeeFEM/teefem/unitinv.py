# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 20:25:52 2012

@author: Jukka

Unit interval

"""

import geom
import numpy as np

def unitinterval(x0=0,x1=1,numberofnodes=5,dense=0):
    '''
    
    Line mesh with linear Seg2 elements

    dense = 0: Unit interval
    x0---o---o---o---x1

    dense > 0: Biased mesh to x1
    x0-----o----o--o-x1

    dense < 0: Biased mesh to x0
    x0-o--o----o-----x1
    
    '''
    
    x = biasedspace(x0,x1,numberofnodes,dense)
    print x
    nodes = [geom.Node(x=xi,y=0,z=0) for xi in x]
    elements = []
    for i in range(numberofnodes-1):
        seg2 = geom.Seg2(nodes=(nodes[i],nodes[i+1]))
        elements.append(seg2)
    mesh = geom.Mesh()
    mesh.elements = elements
    mesh.nodes = nodes
    return mesh

def biasedspace(x0=0,x1=1,numberofnodes=5,dense=1):
    ab0=np.linspace(0,np.log2(2**abs(dense))+1,numberofnodes)
    if np.sign(dense)<>0:
        ab0=2**ab0-1
        ab0=ab0/max(ab0)
        ab0=np.sort(abs(np.sign(dense)-ab0))
        ab0=ab0-min(ab0)
    return (x1-x0)*ab0+x0

def test1(x0=0, x1=1, numberofnodes=5, dense=0):
    mesh = unitinterval(x0,x1,numberofnodes,dense)
    print mesh.nodes
    print mesh.elements
    # => [Node : (0.00,0.00,0.00), Node : (0.25,0.00,0.00), Node : (0.50,0.00,0.00), Node : (0.75,0.00,0.00), Node : (1.00,0.00,0.00)]
    # => [<geom.Seg2 object at 0x046E23D0>, <geom.Seg2 object at 0x046E2490>, <geom.Seg2 object at 0x046E25F0>, <geom.Seg2 object at 0x046E2750>]

if __name__ == '__main__':
    test1(0, 1, 5, -2)
