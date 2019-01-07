# -*- coding: utf-8 -*-

from MoldenData import MoldenData
from Find_coeff_test import findorbs
from Find_coeff_test import findorbsthresh

MD = MoldenData("testfiles/porphyrazine.molden")

#print(MD.basis_labels.index('Fe[9]-4Dzz'))



#print("3Dzz:")
#print(findorbs(MD,"Fe[9]-3Dzz", 10))
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
print(findorbsthresh(MD,"Fe[9]-4Dzz", 0.1))