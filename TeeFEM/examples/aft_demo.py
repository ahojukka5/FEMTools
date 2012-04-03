# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 06:38:31 2012

@author: Jukka Aho

AFT 

"""

from __future__ import division
import numpy as np
import matplotlib.pylab as plt
import time

###############################################################################

# http://www.bryceboe.com/wordpress/wp-content/uploads/2006/10/intersect.py

def ccw(A,B,C):
	return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
	return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

###############################################################################


L = 2.0
r = 1.0

def dist(node1,node2):
    return np.sqrt((node2.x-node1.x)**2 + (node2.y-node1.y)**2)

class Node(object):
    
    count = 0
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.delta = 1
        self.active = True
        self.id = Node.count
#        node.stretch = 1
#        node.nx = 1
#        node.ny = 0 
        Node.count += 1

    def plot_circle(self,r):
        phi = np.linspace(0,2*np.pi)
        x = self.x + r*np.cos(phi)
        y = self.y + r*np.sin(phi)
        plt.plot(x,y,'--b')
        plt.scatter([self.x],[self.y])

node = Node

n1 = node(L-r,0)
n2 = node(L-r/np.sqrt(2),r/np.sqrt(2))
n3 = node(L,r)
n4 = node(L,L)
n5 = node(0.5*L,L)
n6 = node(0,L)
n7 = node(0,0.5*L)
n8 = node(0,0)

class Edge(object):
    count = 0
    def __init__(self, n1, n2):
        self.node1 = n1
        self.node2 = n2
        self.active = True
        self.id = Edge.count
        Edge.count += 1
        
    @property
    def length(self):
        return dist(self.node1, self.node2)
    @property
    def xcenter(self):
        return 0.5*(self.node1.x + self.node2.x)
    @property
    def ycenter(self):
        return 0.5*(self.node1.y + self.node2.y)
    @property
    def angle(self):
        return np.arctan2(self.node2.y-self.node1.y,self.node2.x-self.node1.x)
    @property
    def Kideal(self):
        ''' Returns the Kid of edge '''
        lid = self.length*np.sqrt(3)/2
        aid = self.angle
        xnew = self.xcenter + lid*np.cos(aid+np.pi/2)
        ynew = self.ycenter + lid*np.sin(aid+np.pi/2)
        return Node(xnew,ynew)
    @property
    def delta(self):
        return 0.5*(self.node1.delta + self.node2.delta)
    @property
    def radius(self):
        l = self.length
        return 0.8*min(max(self.delta*l,0.55*l),2.0*l)
    def plot(self):
        plt.plot([self.node1.x,self.node2.x],[self.node1.y,self.node2.y],'r')

edge = Edge

e1 = edge(n1,n2)
e2 = edge(n2,n3)
e3 = edge(n3,n4)
e4 = edge(n4,n5)
e5 = edge(n5,n6)
e6 = edge(n6,n7)
e7 = edge(n7,n8)
e8 = edge(n8,n1)


class Front(object):
    
    def __init__(self, *args):
        self.edges = list(args)
        self.nodes = []
        # Grid-cell parameters
        for edge in self.edges:
            self.nodes.append(edge.node1)
            self.nodes.append(edge.node2)
    
    @property
    def active_nodes(self):
        return [node for node in self.nodes if node.active]

    @property
    def passive_nodes(self):
        return [node for node in self.nodes if not node.active]

    @property
    def active_edges(self):
        return sorted([edge for edge in self.edges if edge.active], key = lambda e: e.length)

    @property
    def passive_edges(self):
        return [edge for edge in self.edges if not edge.active]
    
    def active_nodes_in_radius(self,centernode,radius):
        nodes = []
        for node in self.active_nodes:
            if dist(node,centernode) < radius:
                nodes.append(node)
        return sorted(nodes, key = lambda node: dist(node,centernode))
    
    def intersect(self, edge1, edge2):
        return intersect(edge1.node1,edge1.node2,edge2.node1,edge2.node2)
    
    def plot_active_edges(self):
        for edge in self.active_edges:
            plt.plot([edge.node1.x, edge.node2.x],[edge.node1.y,edge.node2.y],'b--')
            plt.text(edge.xcenter, edge.ycenter, edge.id)
    
    def plot_active_nodes(self):
        x = [node.x for node in self.active_nodes]
        y = [node.y for node in self.active_nodes]
        plt.scatter(x,y)
    
    def find_intersection(self,edge):
        for edg in self.active_edges:
            if edg is edge:
                continue
            if self.intersect(edg, edge):
                print("Intersection between segments: ({e1.node1.x:0.2f},{e1.node1.y:0.2f}) -- ({e1.node2.x:0.2f},{e1.node2.y:0.2f}) and ({e2.node1.x:0.2f},{e2.node1.y:0.2f}) -- ({e2.node2.x:0.2f},{e2.node2.y:0.2f})".format(e1 = edg, e2 = edge))
                return edg
        return None
    
    def step(self):
        t0 = time.clock()
        i = 0
        while True:
            # Find the first active edge
            side = self.active_edges[i]
            # Find the ideal point K_{ideal} for equilateral triangle
            Kideal = side.Kideal
            # Define the radius for a circle 
            Kideal.plot_circle(side.radius)
            # Find all nodes that are in radius
            nodes_in_radius = self.active_nodes_in_radius(Kideal, side.radius)
            print nodes_in_radius

            edge1 = edge(side.node1, Kideal)
            edge2 = edge(side.node2, Kideal)

            if not nodes_in_radius:
                # We are not intersecting any edges, so just add a new edge and make
                # old one passive
                side.active = False
                self.edges.append(edge1)
                self.edges.append(edge2)
                self.nodes.append(Kideal)
                return time.clock()-t0

            intersect = self.find_intersection(edge1) or self.find_intersection(edge2)
            if intersect:# and nodes_in_radius:
#                intersect.plot()
#                side.plot()
                # We are intersecting some edges. Find intersecting edge, and move
                # Kideal to nearest node.
                print("Intersection. Removing edge id {0}: ({n1.x:5.2f},{n1.y:5.2f}) -- ({n2.x:5.2f},{n2.y:5.2f})".format(intersect.id, n1 = intersect.node1, n2 = intersect.node2))
                # Find the "proper" node, it's the nearest one
                Knearest = nodes_in_radius[0] 
                # We must define a passive node and remove it from the front. The passive node is
                # the node that is not the "nearest" node of intersecting line
                passive_node = intersect.node1 if intersect.node1 is not Knearest else intersect.node2
                passive_node.active = False
                # The intersected edge must be passive. Also, we create new edges
                # from side edge to new "not-so-ideal" node just in case that the new
                # edge is not connecting to passive node!
                if side.node1 is not passive_node:
                    self.edges.append( edge(side.node1, Knearest) )
                if side.node2 is not passive_node:
                    self.edges.append( edge(Knearest, side.node2) )
                intersect.active = False
                return time.clock()-t0
            
            # We have failed to create a new edges, so increase i by one and find
            # second shortest edge
            i += 1


# 1
front = Front(e1,e2,e3,e4,e5,e6,e7,e8)
#front.plot_boundary_edges()
#front.plot_boundary_nodes()
itermax = 3
i = 0
while i < itermax:
    i += 1
    plt.figure()
    plt.subplot(111, aspect='equal')
    front.step()
    front.plot_active_edges()
    front.plot_active_nodes()
plt.show()
