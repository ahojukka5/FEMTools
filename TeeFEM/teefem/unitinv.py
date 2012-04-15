# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 20:25:52 2012

@author: Jukka

Unit interval

"""

import geom
import numpy as np

def unitinterval(x0=0,x1=1,numberofnodes=5):
    '''
    
    Unit interval mesh with linear Seg2 elements
    
    x0---o---o---o---x1
    '''
    
    x = np.linspace(x0,x1,numberofnodes)
    nodes = [geom.Node(x=xi,y=0,z=0) for xi in x]
    elements = []
    for i in range(numberofnodes-1):
        seg2 = geom.Seg2(nodes=(nodes[i],nodes[i+1]))
        elements.append(seg2)
    mesh = geom.Mesh()
    mesh.elements = elements
    mesh.nodes = nodes
    return mesh

def test1():
    mesh = unitinterval()
    print mesh.nodes
    print mesh.elements
    # => [Node : (0.00,0.00,0.00), Node : (0.25,0.00,0.00), Node : (0.50,0.00,0.00), Node : (0.75,0.00,0.00), Node : (1.00,0.00,0.00)]
    # => [<geom.Seg2 object at 0x046E23D0>, <geom.Seg2 object at 0x046E2490>, <geom.Seg2 object at 0x046E25F0>, <geom.Seg2 object at 0x046E2750>]

if __name__ == '__main__':
    test1()
