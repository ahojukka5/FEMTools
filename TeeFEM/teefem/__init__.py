# -*- coding: utf-8 -*-

from __future__ import division

import meta


import os
basedir = os.path.dirname(__file__)
datadir = os.path.join(basedir,'data')

#from numpy import *
from common import *
from deprecation import *
from cache import *
from geom import *
from models import *
from materials import *
from elements import *
from boundary_conditions import *
from plate_functions import *
#from scipy.sparse.linalg import spsolve as solve
#
import common
import mesh
#import geom
#import models
#import materials
#import elements
#from scipy.sparse.linalg import spsolve as solve

# "TUI"
#log = common.log
#mesh = geom.mesh
#dkt = models.dkt.dkt

log.info("TeeFEM v{0} loaded".format(meta.__version__))
