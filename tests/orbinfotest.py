#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 01:44:35 2019

@author: mikulas
"""
from OrbInfo import orbinfo
from MoldenData import MoldenData

file = MoldenData("testfiles/H2O.molden")
print(orbinfo(file,2,10))