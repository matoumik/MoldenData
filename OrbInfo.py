#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 00:38:51 2019

@author: mikulas
"""


def orbinfo(moldendatafile, orbnum, num=1):
    orbital = [(x,y) for y,x in sorted(zip(moldendatafile.molecular_orbitals[orbnum], moldendatafile.basis_labels), reverse = True, key = lambda x: abs(x[0]))]
    return orbital[0:num]