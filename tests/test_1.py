import sys
## sys.path.append("D:Ramon_2015/RECERCA/RETOS-2015/Tareas/Proves-FreeCAD-2") # change for your path
#sys.path.append("")
#sys.path.append("/usr/lib/freecad/lib")
import otsun
## import my_materials
import FreeCAD
from FreeCAD import Base
import Part
import time
import numpy as np
## reload(raytrace)

MyProject = 'test_PTC.FCStd'
FreeCAD.openDocument(MyProject)



# ---
# Materials
# ---
otsun.create_simple_volume_material("Glass1", 1.473, 0.015)
otsun.create_reflector_specular_layer("Mir", 0.885, 4.4)
otsun.create_two_layers_material("Mir1", "Mir", "Mir")
otsun.create_absorber_TW_model_layer("Abs1", 0.9427, 0.017, 1.8)

# ---
# Inputs for Total Analysis
# ---

doc = FreeCAD.ActiveDocument
phi_ini = 90.0 + 1.E-9
phi_end = 90.0 + 1.E-4
phi_step = 5.0
theta_ini = 0.0 + 1.E-9
theta_end = 45.0 + 1.E-4
theta_step = 45.0
number_of_rays = 100
aperture_collector_Th = 1845.0 * 10347.0
aperture_collector_PV = 0.0
# for direction of the source two options: Buie model or main_direction 
# direction_distribution = None # default option main_direction
CSR = 0.05
Buie_model = otsun.BuieDistribution(CSR)
direction_distribution = Buie_model
# for integral results three options: ASTMG173-direct (default option), ASTMG173-total, upload data_file_spectrum
data_file_spectrum = 'ASTMG173-direct.txt'
# --------- end

# ---
# Constant inputs for Total Analysis
# ---
show_in_doc = None
polarization_vector = None
light_spectrum = otsun.create_CDF_from_PDF(data_file_spectrum)
# --------- end

power_emitted_by_m2 = otsun.integral_from_data_file(data_file_spectrum)

# objects for scene
sel = doc.Objects
current_scene = otsun.Scene(sel)
results = []
for ph in np.arange(phi_ini, phi_end, phi_step):
    for th in np.arange(theta_ini, theta_end, theta_step):
        main_direction = otsun.polar_to_cartesian(ph, th) * -1.0 # Sun direction vector
        emitting_region = otsun.SunWindow(current_scene, main_direction)
        l_s = otsun.LightSource(current_scene, emitting_region, light_spectrum, 1.0, direction_distribution, polarization_vector)
        exp = otsun.Experiment(current_scene, l_s, number_of_rays, show_in_doc)
        exp.run(show_in_doc)
        t1 = time.time()
        if aperture_collector_Th != 0.0:
            efficiency_from_source_Th = (exp.captured_energy_Th /aperture_collector_Th) / (exp.number_of_rays/exp.light_source.emitting_region.aperture)
        else:
            efficiency_from_source_Th = 0.0
        if aperture_collector_PV != 0.0:
            efficiency_from_source_PV = (exp.captured_energy_PV /aperture_collector_PV) / (exp.number_of_rays/exp.light_source.emitting_region.aperture)
        else:
            efficiency_from_source_PV = 0.0
        results.append((ph, th, efficiency_from_source_Th, efficiency_from_source_PV))

FreeCAD.closeDocument(FreeCAD.ActiveDocument.Name)

def test_1():
    assert 0.9 > results[0][2] > 0.6 and 0.7 > results[1][2] > 0.4 and results[0][3] == 0.0 and results[1][3] == 0.0
