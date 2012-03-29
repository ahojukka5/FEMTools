# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:10:29 2012

Laattafunktioita

@author: Jukka
"""

class PlateCharacteristic(object):
    def __init__(self,**kwds):
        self.thickness = kwds.get('thickness',1)
        self.Tx = kwds.get('Tx', lambda k,e: 0)
        self.Ty = kwds.get('Ty', lambda k,e: 0)
        self.Txy = kwds.get('Txy', lambda k,e: 0)
        self.m = kwds.get('m', lambda k,e: 0)
        self.c = kwds.get('winkler_c', 0)
        
platechar = PlateCharacteristic