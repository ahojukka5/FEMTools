# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 16:45:50 2012

@author: Jukka Aho

######################################
### materials.py - MATERIAL MODELS ###
######################################
"""

import numpy as np

class ElasticIsotropic(object):
    """
    Linear elastic material model
    
    INPUT:
            
    - ``E`` - Modulus of elasticity
    - ``nu`` - Poisson's ratio
        
    OUTPUT:
        
    - Material object
            
    EXAMPLES:
      
    >>> mesh = Mesh(...) # doctest: +SKIP
    >>> mat = Elastic(E = 210e9, nu=0.3)
    >>> mesh.assign_material(mat) # doctest: +SKIP
    
    """
    def __init__(self,*args,**kwds):
        self.type = 'isotropic'
        self.E = kwds.get('E')
        self.nu = kwds.get('nu')
        self.G = self.E/(2*(1+self.nu))

class ElasticOrthotropic(object):
    """
    Linear elastic material model
    
    INPUT:
            
    - ``E`` - Modulus of elasticity
    - ``nu`` - Poisson's ratio
        
    OUTPUT:
        
    - Material object
            
    EXAMPLES:
      
    """
    def __init__(self,*args,**kwds):
        self.type = 'orthotropic'
        self.Ex = kwds.get('Ex')
        self.Ey = kwds.get('Ey')
        self.nux = kwds.get('nux')
        self.nuy = kwds.get('nuy')
        self.G = np.sqrt(self.Ex*self.Ey)/(2*(1+np.sqrt(self.nux*self.nuy)))

elastic = ElasticIsotropic
elastic_orthotropic = ElasticOrthotropic