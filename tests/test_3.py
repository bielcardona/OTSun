import sys
import raytrace
import os
import FreeCAD
from FreeCAD import Base
import Part
import numpy as np

doc = FreeCAD.ActiveDocument


# Materials Solar Cell Stack
raytrace.create_transparent_simple_layer("Trans", 1.0)
raytrace.create_absorber_simple_layer("Abs",1.0)
raytrace.create_two_layers_material("Trans_Abs","Trans","Abs")
file_thin_film = 'Fitxer_OTSun_Exp1a_theta0_90.txt'
file_Perovskite = 'Perovskite_Leguy.txt'
raytrace.create_polarized_thin_film("ThinFilm", file_thin_film, "Vacuum", file_Perovskite)
raytrace.create_PV_material("PV", file_Perovskite)
file_Spiro = 'Spiro_.txt'
raytrace.create_wavelength_volume_material("Spiro", file_Spiro)
file_Ag = 'Ag_Yang.txt'
raytrace.create_metallic_lambertian_layer("Ag", file_Ag)


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



# ---
# Magnitudes used for outputs in Spectral Analysis
# ---
captured_energy_PV = 0.0
captured_energy_Th = 0.0
source_wavelength = []
Th_energy = []
Th_wavelength = []
Th_points_absorber = []
PV_energy = []
PV_wavelength = []
PV_values = []
# --------- end

# objects for scene
sel = doc.Objects
current_scene = raytrace.Scene(sel)

for w in np.arange(wavelength_ini, wavelength_end , wavelength_step):
    light_spectrum = w
    main_direction = raytrace.polar_to_cartesian(phi, theta) * -1.0 # Sun direction vector
    emitting_region = raytrace.SunWindow(current_scene,main_direction)
    l_s = raytrace.LightSource(current_scene, emitting_region, light_spectrum, 1.0, direction_distribution, polarization_vector)
    exp = raytrace.Experiment(current_scene, l_s, number_of_rays, show_in_doc)
    exp.run(show_in_doc)
    print ("%s" % (w)+ '\n')
    Th_energy.append(exp.Th_energy)
    Th_wavelength.append(exp.Th_wavelength)        
    PV_energy.append(exp.PV_energy)
    PV_wavelength.append(exp.PV_wavelength)
    source_wavelength.append(w)       
    if exp.PV_values:
        PV_values.append(exp.PV_values)
    if exp.points_absorber_Th:
        Th_points_absorber.append(exp.points_absorber_Th)
    captured_energy_PV += exp.captured_energy_PV
    captured_energy_Th += exp.captured_energy_Th

# ---
# Output file for wavelengths emitted by the source
# ---
#data_source_wavelength = np.array(np.concatenate(source_wavelength))
data_source_wavelength = np.array(source_wavelength)
data_source_wavelength = data_source_wavelength.T
# --------- end

# ---
# Output source spectrum for calculation and total energy emitted
# ---
source_spectrum = raytrace.spectrum_to_constant_step(data_file_spectrum, wavelength_step, wavelength_ini, wavelength_end)
energy_emitted = np.trapz(source_spectrum[:,1], x = source_spectrum[:,0])
# --------- end


# ---
# Outputs for thermal absorber materials (Th) in Spectral Analysis
# ---
if captured_energy_Th > 1E-9:
    data_Th_points_absorber = np.array(np.concatenate(Th_points_absorber))
    table_Th = raytrace.make_histogram_from_experiment_results(Th_wavelength, Th_energy, wavelength_step, aperture_collector_Th, exp.light_source.emitting_region.aperture)
    spectrum_by_table_Th = source_spectrum[:,1] * table_Th[:,1]		
    power_absorbed_from_source_Th = np.trapz(spectrum_by_table_Th, x = source_spectrum[:,0])
    efficiency_from_source_Th = power_absorbed_from_source_Th / energy_emitted

    # print power_absorbed_from_source_Th * aperture_collector_Th * 1E-6, energy_emitted * exp.light_source.emitting_region.aperture * 1E-6, efficiency_from_source_Th

# --------- end		
	
# ---
# Outputs for photovoltaic materials (PV) in Spectral Analysis
# ---
if captured_energy_PV > 1E-9:

    data_PV_values = np.array(np.concatenate(PV_values))
    table_PV = raytrace.make_histogram_from_experiment_results(PV_wavelength, PV_energy, wavelength_step, aperture_collector_PV, exp.light_source.emitting_region.aperture)
    spectrum_by_table_PV = source_spectrum[:,1] * table_PV[:,1]		
    power_absorbed_from_source_PV = np.trapz(spectrum_by_table_PV, x = source_spectrum[:,0])
    efficiency_from_source_PV = power_absorbed_from_source_PV / energy_emitted
    iqe = internal_quantum_efficiency
    SR = raytrace.spectral_response(table_PV, iqe)
    ph_cu = raytrace.photo_current(SR, source_spectrum)

    # print power_absorbed_from_source_PV * aperture_collector_PV * 1E-6, energy_emitted * exp.light_source.emitting_region.aperture * 1E-6, efficiency_from_source_PV, ph_cu

# --------- end

def test_3():

    assert 0.1 > table_Th[0][1] > 0.0 and 0.1 > table_Th[1][1] > 0.0 and 0.98 > table_PV[0][1] > 0.75 and 0.98 > table_PV[1][1] > 0.75

