# -*- coding: utf-8 -*-

from MoldenData import MoldenData
from pyqint import PyQInt
import numpy as np
from math import sqrt

MD = MoldenData("testfiles/orca_H2O.molden.input", orcacorrect = True)
integrator =  PyQInt()
MD.orca_correct_contractions()

MD.write_molden("testfiles/orca_H2O_fixed.molden.input")


