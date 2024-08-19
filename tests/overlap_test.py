# -*- coding: utf-8 -*-

from MoldenData import MoldenData
from pyqint import PyQInt
import numpy as np
from math import sqrt

MD = MoldenData("testfiles/orca_H2O.molden.input", orcacorrect = True)
#MD_ref = MoldenData("testfiles/h2o_dz.molden")


integrator =  PyQInt()

#print("4S:")
#print(MD.GTO_Overlap("O[1]-2S","O[1]-2S", integrator))
#print(MD_ref.GTO_Overlap("O[1]-2S","O[1]-2S", integrator))


MD.orca_correct_contractions()

MD.calc_GTO_Overlap_matrix()

#MD_ref.calc_GTO_Overlap_matrix()

#print(MD.overlap)
np.savetxt("overlap.txt",MD.overlap,fmt="%10.6f")

v1 = MD.molecular_orbitals[6]
v2 = MD.molecular_orbitals[2]

np.savetxt("v1.txt",v1,fmt="%10.6f")

print(MD.overlap.dot(v2).dot(v1))
#print(v1)

