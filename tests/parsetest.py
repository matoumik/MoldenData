#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 22:24:33 2019

@author: mikulas
"""

from MoldenData import MoldenData

MD = MoldenData()
BG = MoldenData()


MD.MoldenParse("testfiles/H2O.molden", nobeta=False)
print(len(MD.molecular_orbitals))
BG.MoldenParse("testfiles/H2O.molden")
print(len(BG.molecular_orbitals))


#print(MD.atoms)
#print(MD.atom_charges)
#print(len(MD.basis_labels))
#print(len(MD.molecular_orbitals))