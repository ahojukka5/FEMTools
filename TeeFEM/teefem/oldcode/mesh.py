# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:44:29 2012

@author: Jukka Aho
"""

import os
import re
from nets import Node, Seg2, Tria3, Seg3, Tria6, Quad9, Poi1, Quad8, Quad4

import logging
logging.basicConfig(level = logging.DEBUG, filename='mesh.log')

class Mesh(object):
    
    ''' Mesh class. Reads following file formats: Gmsh .msh, Salome .med, Aster .mail '''    
    
    def __init__(self, *args, **kwds):
        ''' Initializes mesh class '''
        logging.debug('Creating mesh')
        self.group_no = {}
        self.group_ma = {}
        self.nodes = {}
        self.nets = {}
        if kwds.has_key('filename'):
            fileformat = os.path.splitext(kwds['filename'])[-1]
            func = {
                '.msh': self.readmshfile, 
#                '.med': self.readmedfile, 
                '.mail': self.readmailfile,
                }
            logging.debug("Reading mesh from file {0}".format(kwds['filename']))
            func[fileformat](*args, **kwds)

    def readmailfile(self, *args, **kwds):
        
        ''' Read Code Aster .mail file format '''

        def titre(*args):
            logging.debug(args)

        def coor_2d(*args):
            nbr = 3
            res = int(len(args)/nbr)
            logging.debug("Found {0} nodes (2D)".format(res))
            for i in xrange(res):
                self.nodes[args[nbr*i]] = Node(x = float(args[nbr*i+1]), y = float(args[nbr*i+2]))

        def coor_3d(*args):
            nbr = 4
            res = int(len(args)/nbr)
            logging.debug("Found {0} nodes (3D)".format(res))
            for i in xrange(res):
                self.nodes[args[nbr*i]] = Node(x = float(args[nbr*i+1]), y = float(args[nbr*i+2]), z = float(args[nbr*i+3]))

        def init_net(**kwds):
            params = kwds['params']
            items = kwds['items']
            nbr = params['nbr']
            res = int(len(items)/nbr)
            logging.debug("Found {0} {1} nets".format(res, kwds['netname']))
            for i in xrange(res):
                netname = items[nbr*i]
                nodenames = items[nbr*i+1:nbr*(i+1)]
                self.nets[netname] = params['netf'](name = netname, nodes = [self.nodes[ni] for ni in nodenames])

        def group_no(*args):
            name = args[0]
            logging.debug("New node group: {0}".format(name))
            self.group_no[name] = set([self.nodes[ni] for ni in args[1:]])

        def group_ma(*args):
            name = args[0]
            logging.debug("New nets group: {0}".format(name))
            self.group_ma[name] = set([self.nets[ni] for ni in args[1:]])
    
        def _F(**kwds):
            return kwds

        functions = {
            'TITRE' : titre,
            'COOR_2D': coor_2d,
            'COOR_3D': coor_3d,
            'GROUP_NO': group_no,
            'GROUP_MA': group_ma,
        }

        functions2 = {
            'POI1' : _F(netf=Poi1, nbr=2),
            'SEG2' : _F(netf=Seg2, nbr=3),
            'SEG3' : _F(netf=Seg3, nbr=4),
            'TRIA3' : _F(netf=Tria3, nbr=4),
            'TRIA6' : _F(netf=Tria6, nbr=7),
            'QUAD4' : _F(netf=Quad4, nbr=5),
            'QUAD8' : _F(netf=Quad8, nbr=9),
            'QUAD9' : _F(netf=Quad9, nbr=10),
        }
        with open(kwds['filename'], 'r') as meshfile: data = meshfile.read()
        for section in re.split('FINSF', data):
            items = re.split('\s+', section)
            items = [i.strip('% ') for i in items]
            items = [i for i in items if len(i) != 0]
            if items[0] in functions:
                functions[items[0]](*items[1:])
            if items[0] in functions2:
                init_net(items = items[1:], netname = items[0], params = functions2[items[0]])

    def readmshfile(self, *args, **kwds):
        
        ''' Read Gmsh .msh file format '''
        
        mapping = {
            1: Seg2, 
            2: Tria3, 
            3: Quad4,
            8: Seg3,
            9: Tria6,
            10: Quad9,
            15: Poi1,
            16: Quad8,
            }
        regions = {}
        with open(kwds['filename']) as mshfile:
            section = None
            for line in mshfile:
                line = line.strip()
                if line.startswith('$'):
                    if line.startswith('$End'): section = None
                    else: section = line
                    continue
                columns = line.split()
                if len(columns)==1: continue # Emme tarvitse lukumäärää                
                if section == "$MeshFormat": self.meshformat = line
                if section == "$PhysicalNames":
                    name = columns[2].strip('"')
                    self.group_ma[name] = set()
                    regions[int(columns[1])] = self.group_ma[name]
                    logging.debug("New group_ma: {0}".format(name))
                if section == "$Nodes":
                    self.nodes['N%d'%(int(columns[0]))] = Node(x = float(columns[1]),y = float(columns[2]),z = float(columns[3]))
                if section == "$Elements":
                    columns = [int(ci) for ci in columns]
                    eid,etypeid,enumofparams,eregion,egentity = columns[0:5]
#                    eparams = columns[4:enumofparams-2]
                    enodes = columns[enumofparams+3:]
#                    self.elements[eid] = mapping[etypeid](self.regions[eregion],egentity,eparams,[self.nodes[ni] for ni in enodes])
                    try:
                        net = mapping[etypeid](name = 'M%d'%(int(eid)), nodes = [self.nodes['N%d'%(ni)] for ni in enodes])
                    except KeyError:
                        print("Unable to create net type {0}".format(etypeid))
                        print columns
                    self.nets['M%d'%(int(eid))] = net
                    regions[eregion].add(net)
        self.create_node_groups()

#    def readmedfile(self, *args, **kwds): # HARD LINK
#
#       ''' Read Salome .med file format '''
#
#        print("I'm broken")
#        return 0
#        db = tables.openFile(kwds['filename'])
#        NOE = db.getNode('/ENS_MAA/geom/NOE')
#        MAI = db.getNode('/ENS_MAA/geom/MAI')
#        FAS = db.getNode('/ENS_MAA/geom/FAS')
#        # GROUP_MA
#        for group_ma in FAS.ELEME._f_iterNodes():
#            name = ''.join(group_ma.GRO.NOM[:]).strip()
#            idx = group_ma._v_attrs.NUM
#            GMA = GROUP_MA(name, idx)
#            console_print("New group_ma: {0} with index {1}".format(name,idx))
#            self.group_ma_idx[idx] = GMA
#            self.group_ma_name[name] = GMA
#        # NODES
#        COO = NOE.COO[:]
#        FAM = NOE.FAM[:]
#        NBR = NOE.COO._v_attrs.NBR
#        print("Number of nodes: %d"%(NBR))
#        for idx in xrange(NBR):
#            fam = FAM[idx]
#            coo = COO[idx::NBR]
#            node = Node(coo)
#            self.nodes[idx+1] = node
#        # NETS
#        mapping = {
#            "PO1": POI1, 
#            "SE2": SEG2,
#            "SE3": SEG3,
#            "TR3": TRIA3,
#            "TR6": TRIA6,
#            }
#        for (typ,arrs) in MAI._v_children.iteritems():
#            NBR = int(arrs.NOD._v_attrs.NBR)
#            console_print("Number of {type} nets: {nbr:d}".format(type=typ, nbr=NBR))
#            ntyp = mapping[typ]
#            NOD = arrs.NOD[:]
#            FAM = arrs.FAM[:]
#            for idx in xrange(NBR):
#                fam = FAM[idx]
#                nod = NOD[idx::NBR]
##                net = ntyp(nodes = [self.nodes[ni] for ni in nod][::-1], group_ma = self.group_ma_idx[fam])
#                net = ntyp(nodes = [self.nodes[ni] for ni in nod], group_ma = self.group_ma_idx[fam])
#                self.nets.add(net)
        
    def create_node_groups(self):
        ''' Create node groups from element groups '''
        for key in self.group_ma.iterkeys():
            self.group_no[key] = set()
            logging.debug("New node group: {0}".format(key))
            for element in self.group_ma[key]:
                for node in element.nodes:
                    self.group_no[key].add(node)

    def assign_material(self, mat):
        for n in self.nets.itervalues():
            n.material = mat
            