from .math import random_congruential
import time
import numpy as np


class Experiment:
    """
    Sets up and runs and experiment in a given scene with a given light source.
    If show_in_doc is given, the emitting region is drawn in the FreeCAD active document.
    If show_in_doc is given, the rays could be drawn in the FreeCAD active document (using the run function).
    """

    def __init__(self, scene, light_source, number_of_rays, show_in_doc=None):
        self.scene = scene
        self.light_source = light_source
        if show_in_doc:
            self.light_source.emitting_region.add_to_document(show_in_doc)
        self.number_of_rays = number_of_rays
        self.wavelengths = []
        self.captured_energy_Th = 0
        self.captured_energy_PV = 0
        self.Th_energy = []
        self.Th_wavelength = []
        self.PV_energy = []
        self.PV_wavelength = []
        self.PV_values = []
        self.points_absorber_Th = []
        random_congruential(time.time()) # TODO: change location

    def run(self, show_in_doc=None):
        for _ in np.arange(0,self.number_of_rays,1):
            ray = self.light_source.emit_ray()
            ray.run()
            self.wavelengths.append(ray.wavelength)
            if show_in_doc:
                ray.add_to_document(show_in_doc)
            if ray.got_absorbed:
                self.captured_energy_Th += ray.energy
                self.Th_energy.append(ray.energy)
                self.Th_wavelength.append(ray.wavelength)
                self.points_absorber_Th.append((ray.energy,
                                             ray.points[-1].x, ray.points[-1].y, ray.points[-1].z,
                                             ray.points[-2].x, ray.points[-2].y, ray.points[-2].z,
                                             ray.normals[-1].x, ray.normals[-1].y, ray.normals[-1].z))
            else:
                self.Th_energy.append(0.0)
                self.Th_wavelength.append(ray.wavelength)
            if ray.in_PV:
                PV_energy_absorbed = np.sum(ray.PV_absorbed)
                self.captured_energy_PV += PV_energy_absorbed
                self.PV_energy.append(PV_energy_absorbed)
                self.PV_wavelength.append(ray.wavelength)
                length = len(ray.PV_values)
                if length > 0:
                    for z in np.arange(0,length,1):
                        self.PV_values.append(ray.PV_values[z])
            else:
                self.PV_energy.append(0.0)
                self.PV_wavelength.append(ray.wavelength)
