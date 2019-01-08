# -*- coding: utf-8 -*-

from MoldenData import MoldenData
from Find_coeff_test import findorbs
from Find_coeff_test import findorbsthresh

MD = MoldenData("FePz_q.molden")

d = ("zz", "xy", "yz", "xz", "xx-yy")
O3D = ["Fe[9]-3D" + o for o in d]
O4D = ["Fe[9]-4D" + o for o in d]
O5D = ["Fe[9]-5D" + o for o in d]

thresh = 0.05
print("\n\n\n\n\n")
print("Normovane")

print("3D")
print(findorbs(MD ,O3D , 30))
print("4D")
print(findorbs(MD ,O4D , 30))
print("5D")
print(findorbs(MD ,O5D , 30))


print("3Dzz:")
print(findorbs(MD,"Fe[9]-3Dzz", 10))
print("4Dzz:")
print(findorbs(MD,"Fe[9]-4Dzz", 10))
print("5Dzz:")
print(findorbs(MD,"Fe[9]-5Dzz", 10))

print("3Dxy:")
print(findorbs(MD,"Fe[9]-3Dxy", 10))
print("4Dxy:")
print(findorbs(MD,"Fe[9]-4Dxy", 10))
print("5Dxy:")
print(findorbs(MD,"Fe[9]-5Dxy", 10))

print("3Dyz:")
print(findorbs(MD,"Fe[9]-3Dyz", 10))
print("4Dyz:")
print(findorbs(MD,"Fe[9]-4Dyz", 10))
print("5Dyz:")
print(findorbs(MD,"Fe[9]-5Dyz", 10))

print("3Dxz:")
print(findorbs(MD,"Fe[9]-3Dxz", 10))
print("4Dxz:")
print(findorbs(MD,"Fe[9]-4Dxz", 10))
print("5Dxz:")
print(findorbs(MD,"Fe[9]-5Dxz", 10))


print("3Dxx-yy:")
print(findorbs(MD,"Fe[9]-3Dxx-yy", 10))
print("4Dxx-yy:")
print(findorbs(MD,"Fe[9]-4Dxx-yy", 10))
print("5Dxx-yy:")
print(findorbs(MD,"Fe[9]-5Dxx-yy", 10))

print("4S:")
print(findorbs(MD,"Fe[9]-4S", 10))

print("4Px:")
print(findorbs(MD,"Fe[9]-4Px", 10))
print("4Py:")
print(findorbs(MD,"Fe[9]-4Py", 10))
print("4Pz:")
print(findorbs(MD,"Fe[9]-4Pz", 10))


print("______________________________________________________\n\n")
print("Nenormovane")
print("thresh:   " + str(thresh))
print("\n")

for oname in O4D:
  print(oname)
  print(findorbsthresh(MD,oname,thresh))
