# -*- coding: utf-8 -*-
"""

"""


def findorbs(moldendatafile, orbs, num=1):
    orbindices=list()
    orbcoeffpart=list()
    orbcoeffs=moldendatafile.molecular_orbitals
    orblabels=moldendatafile.basis_labels

    if isinstance(orbs, str):
        try:
            orbindices.append(orblabels.index(orbs))
        except ValueError:
            print("Error: Orb not found")    
    else:
        for o in orbs:
            try:
                orbindices.append(orblabels.index(o))
            except ValueError:
                print("Error: Orb not found")

    
        
        
    for orb in orbcoeffs:
        if len(orb)>0:
            norm=sum( i*i for i in orb)
            coeff=sum(orb[i]*orb[i] for i in orbindices)
            if norm>0:
                orbcoeffpart.append(coeff/norm)
            else: orbcoeffpart.append(-1.)
    
    largest = sorted(range(len(orbcoeffpart)),key=lambda x: orbcoeffpart[x], reverse = True)
    return largest[0:num]

def findorbsthresh(moldendatafile, orbs, thresh):
    orbindices=list()
    orbreturn=list()
    orbcoeffsreturn=list()
    orbcoeffs=moldendatafile.molecular_orbitals
    orblabels=moldendatafile.basis_labels

    if isinstance(orbs, str):
        try:
            orbindices.append(orblabels.index(orbs))
        except ValueError:
            print("Error: Orb not found")    
    else:
        for o in orbs:
            try:
                orbindices.append(orblabels.index(o))
            except ValueError:
                print("Error: Orb not found")
    
    for orb in orbcoeffs:
        if len(orb)>0:
            coeff=sum(orb[i]*orb[i] for i in orbindices)
            if coeff>thresh:
                orbreturn.append(orbcoeffs.index(orb))
                orbcoeffsreturn.append(coeff)
    
    orbreturn = sorted(orbreturn,key=lambda x: orbcoeffs[x], reverse = True)
    return orbreturn