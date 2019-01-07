# -*- coding: utf-8 -*-
"""
This is a class for parsing data from a single molden file.
"""

def todo():
    print("Not yet implemented")    

class MoldenData:
    
    KEYWORDS = ()
    
    def __init__(self, filename=None, nobeta=True):
        
        self.atoms = list()
        self.atom_charges = list()
        self.atom_positions = list()
        self.basis_labels = list()
        self.molecular_orbitals = list()
        self.orbital_info = list()
        
        self.units = "Unknown"
        self.atom_number = 0
        
        
        if filename != None:
            self.MoldenParse(filename,nobeta)
    
    def MoldenParse(self, filename, nobeta=True):
        moldenfile = open(filename, 'r')
        print("Opened file:  "+filename)
        block = ""
        _atom_orb_num = 0 
        _orb_coeff = list()
        _orb_info = dict()
        for line in moldenfile:
            fline = line.split()
            if len(fline)>0:
                if fline[0][0] == '[':
                    block = "".join(fline).upper()
                    block = self._react(block)
                elif block == "[ATOMS]":
                    self.atoms.append(fline[0])
                    self.atom_charges.append(int(fline[2]))
                    self.atom_positions.append(tuple([float(i) for i in fline[2:5]]))
                
                elif block == "[GTO]":
                    if len(fline) > 1 and fline[1] == "0":
                        _s_num = 0
                        _p_num = 1
                        _d_num = 2
                        _f_num = 3
                        _g_num = 4
                        _atom_orb_num +=1;
                    elif fline[0].isalpha():
                        l = fline[0].upper()
                        if l == "S":
                            _s_num += 1
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}S'.format(_atom_orb_num,_s_num))
                        elif l == "P":
                            _p_num += 1
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Px'.format(_atom_orb_num,_p_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Py'.format(_atom_orb_num,_p_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Pz'.format(_atom_orb_num,_p_num))
                        elif l == "D":
                            _d_num += 1
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dzz'.format(_atom_orb_num,_d_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxz'.format(_atom_orb_num,_d_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dyz'.format(_atom_orb_num,_d_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxy'.format(_atom_orb_num,_d_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxx-yy'.format(_atom_orb_num,_d_num))
                        elif l == "F":
                            _f_num += 1
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fzzz'.format(_atom_orb_num,_f_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fxzz'.format(_atom_orb_num,_f_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fyzz'.format(_atom_orb_num,_f_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fxyz'.format(_atom_orb_num,_f_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Fzxx-zyy'.format(_atom_orb_num,_f_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dxxx-3xyy'.format(_atom_orb_num,_f_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}Dyyy-3yxx'.format(_atom_orb_num,_f_num))
                        elif l == "G": 
                            _g_num += 1
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G1'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G2'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G3'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G4'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G5'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G6'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G7'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G8'.format(_atom_orb_num,_g_num))
                            self.basis_labels.append(self.atoms[_atom_orb_num-1]+'[{}]-{}G9'.format(_atom_orb_num,_g_num))
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
        elif keyword == "[ATOMS]ANGS":
            self.units = "ANGS"
            keyword = "[ATOMS]"


        return keyword