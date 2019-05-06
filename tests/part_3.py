import sys
import otsun
import os
import FreeCAD
from FreeCAD import Base
import Part
import numpy as np

MyProject = 'Perovskite_Stack_200nm.FCStd'
FreeCAD.openDocument(MyProject)


doc = FreeCAD.ActiveDocument


# Materials Solar Cell Stack
otsun.TransparentSimpleLayer("Trans", 1.0)
otsun.AbsorberSimpleLayer("Abs", 1.0)
otsun.TwoLayerMaterial("Trans_Abs", "Trans", "Abs")
file_thin_film = 'Fitxer_OTSun_Exp1a_theta0_90.txt'
file_Perovskite = 'Perovskite_Leguy.txt'
otsun.PolarizedThinFilm("ThinFilm", file_thin_film, "Vacuum", file_Perovskite)
otsun.PVMaterial("PV", file_Perovskite)
file_Spiro = 'Spiro_.txt'
otsun.WavelengthVolumeMaterial("Spiro", file_Spiro)
file_Ag = 'Ag_Yang.txt'
otsun.MetallicLambertianLayer("Ag", file_Ag)


# ---
# Constant inputs for Spectral Analysis
# ---
polarization_vector = None
show_in_doc = None
# show_in_doc = doc
# --------- end

# ---
# Inputs for Spectral Analysis
# ---
phi = 0.0 # default value zero
phi = phi + 1.E-9
theta = 0.0 # default value zero
theta = theta + 1.E-9
wavelength_ini = 500.0 # default value 280.0
wavelength_end = 502.0 # default value 4000.0
wavelength_end = wavelength_end + 1E-4
wavelength_step = 2.0 # default value 10.0
number_of_rays = 100 # number of rays per wavelength # default value 1000
aperture_collector_Th = 1000. * 1000. * 1.0 # default value zero
aperture_collector_PV = 1000. * 1000. * 1.0 # default value zero
# for direction of the source two options: Buie model or main_direction 
direction_distribution = None # default option main_direction
# CSR = 0.05
# Buie_model = raytrace.BuieDistribution(CSR)
# direction_distribution = Buie_model
# for the internal quantum efficiency two options: constant value =< 1.0, or data file 
internal_quantum_efficiency = 1.0 # default option equal to 1.0
# internal_quantum_efficiency = 'data.txt'  
data_file_spectrum = 'ASTMG173-direct.txt'
# --------- end


FreeCAD.closeDocument(FreeCAD.ActiveDocument.Name)

