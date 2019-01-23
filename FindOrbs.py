#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 10:52:44 2019

@author: mikulas
"""

def findorbs_coeff(moldendatafile, orbs, num=1):
    """
    Searches the set of molekular orbitals and returns the specified number 
    of orbitals with the highest relative contribution of the basis set function.
    Sorted by the value of 
    (sum_of_squares_of_input_coefficients/sum_of_squares_of_all_coefficients).
    Accepts eiher a string, or a list of strings.
    """
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

def findorbs_max(moldendatafile, orbs):
    """
    Finds all molecular orbitals, in which the coefficient of the selected basis set function 
    is higher than any other basis function coefficient. Accepts eiher a string, or a list of strings.
    In the case of a list, returns all molecular orbitals,
    for which any of the input basis functions is dominant.
    """
    orbindices=list()
    orbsout=list()
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
        if orb.index(max(orb)) in orbindices:
            orbsout.append(orbcoeffs.index(orb))
            
    return orbsout

    