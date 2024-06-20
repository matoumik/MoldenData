# -*- coding: utf-8 -*-
"""
This is a class for parsing data from a single molden file.
"""
from pyqint import PyQInt, cgf
from copy import deepcopy
from math import sqrt,pi
import numpy as np


def todo():
    print("Not yet implemented")    

class MoldenData:
    
    KEYWORDS = ()
    
    def __init__(self, filename=None, nobeta=True, orcacorrect = False):
        
        self.atoms = list()
        self.atom_charges = list()
        self.atom_positions = list()
        self.basis_labels = list()
        self.molecular_orbitals = list()
        self.orbital_info = list()
        self.gto_coeffs = dict()
        
        self.units = "Unknown"
        self.atom_number = 0
        self.overlap = False
        
        
        if filename != None:
            self.MoldenParse(filename,nobeta)
        
        if orcacorrect:
            self.orca_correct_contractions()
    
    def MoldenParse(self, filename, nobeta=True):
        moldenfile = open(filename, 'r')
        print("Opened file:  "+filename)
        block = ""
        _atom_orb_num = 0 
        _orb_coeff = list()
        _orb_info = dict()
        _labels = list()
        for line in moldenfile:
            fline = line.split()
            if len(fline)>0:
                if fline[0][0] == '[':
                    block = "".join(fline).upper()
                    block = self._react(block)
                elif block == "[ATOMS]":
                    self.atoms.append(fline[0])
                    self.atom_charges.append(int(fline[2]))
                    self.atom_positions.append(tuple([float(i) for i in fline[3:6]]))
                    #print(fline)
                
                elif block == "[GTO]":
                    if len(fline) > 1 and fline[1] == "0":
                        _s_num = 0
                        _p_num = 1
                        _d_num = 2
                        _f_num = 3
                        _g_num = 4
                        _atom_orb_num +=1
                    elif fline[0].isalpha():
                        
                        for i in _labels:
                            self.gto_coeffs[i]=deepcopy(_GTOs)
                        
                        _labels = list()
                        _GTOs = list()
                        _lcoeff = 1
                        
                        l = fline[0].upper()
                        if l == "S":
                            _s_num += 1
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}S'.format(_atom_orb_num,_s_num))
                            self.basis_labels += _labels
                        elif l == "P":
                            _p_num += 1
                            _lcoeff = 3
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Px'.format(_atom_orb_num,_p_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Py'.format(_atom_orb_num,_p_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Pz'.format(_atom_orb_num,_p_num))
                            self.basis_labels += _labels
                        elif l == "D":
                            _d_num += 1
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dzz'.format(_atom_orb_num,_d_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxz'.format(_atom_orb_num,_d_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dyz'.format(_atom_orb_num,_d_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxx-yy'.format(_atom_orb_num,_d_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxy'.format(_atom_orb_num,_d_num))
                            self.basis_labels += _labels
                        elif l == "F":
                            _f_num += 1
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fzzz'.format(_atom_orb_num,_f_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fxzz'.format(_atom_orb_num,_f_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fyzz'.format(_atom_orb_num,_f_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fxyz'.format(_atom_orb_num,_f_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fzxx-zyy'.format(_atom_orb_num,_f_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fxxx-3xyy'.format(_atom_orb_num,_f_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fyyy-3yxx'.format(_atom_orb_num,_f_num))
                            self.basis_labels += _labels
                        elif l == "G": 
                            _g_num += 1
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G1'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G2'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G3'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G4'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G5'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G6'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G7'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G8'.format(_atom_orb_num,_g_num))
                            _labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G9'.format(_atom_orb_num,_g_num))
                            self.basis_labels += _labels
                    else:
                        _gto = [float(i.replace('D', 'E')) for i in fline]
                        _GTOs.append(_gto)
                elif block == "[MO]":
                    if fline[0].isdigit():
                        _orb_coeff.append(float(fline[1]))
                        
                        if len(_orb_info)>0:
                            self.orbital_info.append(_orb_info)
                            _orb_info = dict()
                    else:
                        if len(_orb_coeff)>0:
                            if len(_orb_coeff)==len(self.basis_labels):
                                self.molecular_orbitals.append(_orb_coeff)
                                _orb_coeff = list()
                            else:
                                print("Error: Number of orbital coefficients did not match basis size!")
                                _orb_coeff.clear()
                        
                        _templine=" ".join(fline) #Protože někdy z rovnítkem není mezera... :(
                        _templine = _templine.replace("="," ")
                        fline = _templine.split()
                        
                        
                        if fline[0].upper() == "SYM":
                            pass
                        elif fline[0].upper() == "ENE":
                            _orb_info[0] = float(fline[1])
                        elif fline[0].upper() == "SPIN":
                            _orb_info[1] = fline[1].upper()
                        elif fline[0].upper() == "OCCUP":
                            _orb_info[2] = float(fline[1])
                            
                            
        if len(_orb_coeff)>0: #Protože konec souboru nemám jinak moc šanci poznat...
            if len(_orb_coeff)==len(self.basis_labels):
                self.molecular_orbitals.append(_orb_coeff)
                _orb_coeff = list()
            else:
                print("Error: Number of orbital coefficients did not match basis size!")
                _orb_coeff.clear()

        if len(_labels)>0:
            for i in _labels:
                self.gto_coeffs[i]=deepcopy(_GTOs)
                
        print("Finished reading MOs...")
                
        moldenfile.close()
        print("Closed file:  " + filename)
        
        if nobeta:
            self.molecular_orbitals = [orb for orb in self.molecular_orbitals if self.orbital_info[self.molecular_orbitals.index(orb)][1]!="BETA"]
            print("Removed beta orbitals...") 
        
        self.atom_number = len(self.atoms)
    
        
    def SanityCheck(self):
        pass
    
    def _react(self, keyword):
        #print(keyword)
        if keyword == "[ATOMS]AU":
            self.units = "AU"
            keyword = "[ATOMS]"
        elif keyword == "[ATOMS](AU)":
            print("PSI4 file detected")
            self.units = "AU"
            keyword = "[ATOMS]"
        elif keyword == "[ATOMS]ANGS":
            self.units = "ANGS"
            keyword = "[ATOMS]"
        return keyword
        
        
    def GTO_Overlap(self, lab_i, lab_j, integrator):
        atom1 = int(lab_i.split("[")[1].split("]")[0])-1
        atom2 = int(lab_j.split("[")[1].split("]")[0])-1

        if(self.units == "AU"):
            scale = 1.0
        elif(self.units == "ANGS"):
            scale = 1.8897259886
        else:
            print("Error: invalid units for geometry")

        pos1 = [ i*scale for i in self.atom_positions[atom1]]
        pos2 = [ i*scale for i in self.atom_positions[atom2]]

        #print(self.atom_positions)

        #print(pos1)
        #print(pos2)
    	
        cgf1 = cgf(pos1)
        cgf2 = cgf(pos2)
    	
        #print(lab_i," ",lab_j)

        for lab, cgfi in [(lab_i,cgf1),(lab_j,cgf2)]:
            for l in self.gto_coeffs[lab]:
                if "Px" in lab:
    	            cgfi.add_gto(l[1], l[0], 1,0,0)
                elif "Py" in lab:
    	            cgfi.add_gto(l[1], l[0], 0,1,0)
                elif "Pz" in lab:
    	            cgfi.add_gto(l[1], l[0], 0,0,1)
                elif "Dxy" in lab:
                    cgfi.add_gto(l[1], l[0], 1,1,0)
                elif "Dxz" in lab:
                    cgfi.add_gto(l[1], l[0], 1,0,1)
                elif "Dyz" in lab:
                    cgfi.add_gto(l[1], l[0], 0,1,1)
                elif "Dxx-yy" in lab:  
                    cgfi.add_gto(l[1]/sqrt(4/3), l[0], 2,0,0)
                    cgfi.add_gto(-l[1]/sqrt(4/3), l[0], 0,2,0)
                  
                    
                elif "Dzz" in lab:
                    cgfi.add_gto(  l[1], l[0], 0,0,2)
                    cgfi.add_gto( -l[1]/2, l[0], 2,0,0)
                    cgfi.add_gto( -l[1]/2, l[0], 0,2,0)
                    
                elif "Fzzz" in lab:
                    cgfi.add_gto(2*l[1]/sqrt(6), l[0], 0,0,3)
                    cgfi.add_gto( -l[1]/sqrt(6), l[0], 3,0,0)
                    cgfi.add_gto( -l[1]/sqrt(6), l[0], 0,3,0)
                elif "Fxzz" in lab:
                    cgfi.add_gto(4*l[1]/sqrt(18), l[0], 1,0,2)
                    cgfi.add_gto( -l[1]/sqrt(18), l[0], 3,0,0)
                    cgfi.add_gto( -l[1]/sqrt(18), l[0], 1,2,0)
                
                else:
                    #print("Sorb used")
                    cgfi.add_gto(l[1], l[0], 0,0,0)
        return integrator.overlap(cgf1,cgf2)

    	        
   
    def calc_GTO_Overlap_matrix(self):
        print("Calculating overlap matrix:")

        dim = len(self.basis_labels)
        
        integrator = PyQInt()

        self.normalize_CGFs(integrator)

        self.overlap = np.zeros((dim,dim))
        
        for i in range(dim):
            for j in range(dim):
                self.overlap[i][j]=self.GTO_Overlap(self.basis_labels[i],self.basis_labels[j],integrator)

        print("Overlap matrix done")
    
    
    def normalize_CGFs(self, integrator):
        
        for orb in self.basis_labels:
            _norm = self.GTO_Overlap(orb,orb,integrator)
            if(abs(_norm-1)>1e-6):
                print("Renormalizing orbital ", orb, " (norm = " , _norm, ")" )
                for l in self.gto_coeffs[orb]:
                    l[1] = l[1]/sqrt(_norm)

    def orca_correct_contractions(self):
        print("Correcting contraction coefficients from orca")
        for lab,_cgf in self.gto_coeffs.items():
            for _gto in _cgf:
                if ("Px" in lab) or "Py" in lab or "Pz" in lab:
                    _gto[1] = _gto[1]/(_gto[0]*2/pi)**(3/4)/sqrt(_gto[0]*4)                
                elif "Dxy" in lab or "Dxz" in lab or "Dyz" in lab or "Dxx-yy" in lab or "Dzz" in lab:
                    _gto[1] = _gto[1]/(_gto[0]*2/pi)**(3/4)/sqrt(_gto[0]**2*16)
                elif "F" in lab.split("]")[1]:
                    _gto[1] = _gto[1]/(_gto[0]*2/pi)**(3/4)/sqrt(_gto[0]**3*64)
                elif "G" in lab.split("]")[1]:
                    _gto[1] = _gto[1]/(_gto[0]*2/pi)**(3/4)/sqrt(_gto[0]**4*256/3)
                else:
                    _gto[1] = _gto[1]/(_gto[0]*2/pi)**(3/4)
                

    def write_molden(self,filename):
        print("Writing molden file: " + filename)
        file = open(filename,"w")        
        _head = """[Molden Format]
[Title]
 Molden file created by MoldenData script

[Atoms] """ + self.units + "\n"

        file.write(_head)
        for i in range(len(self.atoms)):
            file.write("{0:<2} {1:3d} {2:4d}   {3:20.10f} {4:20.10f} {5:20.10f}\n".format(self.atoms[i],i+1,self.atom_charges[i],self.atom_positions[i][0],self.atom_positions[i][1],self.atom_positions[i][2]))

        file.write("[GTO]\n")
        
        for i in range(len(self.atoms)):
            file.write( "   " + str(i+1) + " 0\n" )
            _num = dict()
            _num["s"]=0
            _num["p"]=0
            _num["d"]=0
            _num["f"]=0
            _num["g"]=0

            for lab,_cgf in self.gto_coeffs.items():
                if ("["+str(i+1)+"]") in lab:
                    _olab = lab.split("]-")[1]
                    _type = "x"
                    if("S" in _olab): 
                        _type = "s"
                    elif("Px" in _olab): 
                        _type = "p"
                    elif("Dzz" in _olab): 
                        _type = "d"
                    elif("Fzzz" in _olab):
                        _type = "f"
                    elif("G1" in _olab):
                        _type = "g"
                    
                    if(_type != "x"):
                        _num[_type]+=1
                        file.write(_type+"  "+str(len(_cgf))+"  "+"1.00"+"\n")
                        for _gto in _cgf:
                            file.write("{0:20.10f} {1:20.10f}\n".format(_gto[0] ,_gto[1]))
            file.write("\n")


        _txt = "[5D]\n[7F]\n[9G]\n[MO]\n"
        file.write(_txt)

        for i in range(len(self.orbital_info)):
            file.write("Ene= " + str(self.orbital_info[i][0]) + "\nSpin= " + self.orbital_info[i][1] + "\nOccup= " + str(self.orbital_info[i][2]) + "\n")

            j = 1
            for _c in self.molecular_orbitals[i]:
                file.write("{0:3d} {1:20.10f}\n".format(j,_c))
                j+=1


      



    	
    	
    
    	
    	
    
