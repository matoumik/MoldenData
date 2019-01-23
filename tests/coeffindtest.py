# -*- coding: utf-8 -*-

from MoldenData import MoldenData
from FindOrbs import findorbs_coeff
from FindOrbs import findorbs_max

MD = MoldenData("testfiles/porphyrazine.molden")

#print(MD.basis_labels.index('Fe[9]-4Dzz'))



#print("3Dzz:")
print(findorbs_coeff(MD,"Fe[9]-3Dzz", 10))
#print("4Dzz:")
#print(findorbs(MD,"Fe[9]-4Dzz", 10))
#print("5Dzz:")
#print(findorbs(MD,"Fe[9]-5Dzz", 10))
#
#print("3Dxy:")
#print(findorbs(MD,"Fe[9]-3Dxy", 10))
#print("4Dxy:")
#print(findorbs(MD,"Fe[9]-4Dxy", 10))
#print("5Dxy:")
#print(findorbs(MD,"Fe[9]-5Dxy", 10))

print("4Dzz:")
print(findorbs_max(MD,"Fe[9]-4S"))