# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 11:10:29 2012

Laattafunktioita

@author: Jukka
"""

class PlateCharacteristic(object):
    def __init__(self,**kwds):
        self.thickness = kwds.get('thickness',1)

platechar = PlateCharacteristic