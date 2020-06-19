"""Module otsun.experiments for creating experiments.

The module defines a class `Experiment` that deals with the setup and
run of experiments.
"""

import numpy as np

class Experiment:
    """
    Sets up and runs and experiment in a given scene with a given light source.

    An Experiment is formed by a scene and a ligh_source. Rays are emitted from light_source
    and interfere with the scene until they are absorbed or scape from the scene. If they are
    absorbed, the energies are computed and stored.

    Parameters
    ----------
    scene : otsun.Scene
        Scene of the experiment
    light_source : otsun.LightSource
        Source of the rays to emit
    number_of_rays : int
        Number of rays to emit in the experiment
    show_in_doc : App.Document
        FreeCAD document where to plot the rays, or None if plotting is not desired

    Attributes
    ----------
    wavelengths: List of float
        List of wavelengths of emitted rays
    captured_energy_Th : float
        Total energy of rays that got absorbed
    captured_energy_PV : float
        Total energy of rays that fell in a PV
    Th_energy :
    Th_wavelength : List of float
        List of wavelengths of rays that got absorbed
    PV_energy : List of floats
        List of energies of rays that fell in a PV
    PV_wavelength : List of float
        List of wavelengths of rays that fell in a PV
    PV_values : List of tuples of floats
        List of PV_values of all emitted rays that fell in a PV
    points_absorber_Th : List of tuple of floats
        List with data (energy, location,...) for each ray that got absorbed
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

    def run(self, show_in_doc=None):
        """Runs the experiment and plots the rays in the document specified (if any)"""
        for _ in np.arange(0, self.number_of_rays, 1):
            ray = self.light_source.emit_ray()
            ray.run()
            self.wavelengths.append(ray.wavelength)
            if show_in_doc:
                ray.add_to_document(show_in_doc)
            if ray.Th_absorbed:
                self.captured_energy_Th += ray.energy
                self.Th_energy.append(ray.energy)
                self.Th_wavelength.append(ray.wavelength)
                self.points_absorber_Th.append((ray.energy,
                                                ray.points[-1].x, ray.points[-1].y, ray.points[-1].z,
                                                ray.points[-2].x, ray.points[-2].y, ray.points[-2].z,
                                                ray.last_normal.x, ray.last_normal.y, ray.last_normal.z))
            else:
                self.Th_energy.append(0.0)
                # TODO: Review... ray.wavelength always added to Th_wavelength
                # Hence always Th_wavelength == wavelenghts
                self.Th_wavelength.append(ray.wavelength)
            if ray.PV_absorbed:
                PV_energy_absorbed = np.sum(ray.PV_absorbed)
                self.captured_energy_PV += PV_energy_absorbed
                self.PV_energy.append(PV_energy_absorbed)
                self.PV_wavelength.append(ray.wavelength)
                length = len(ray.PV_values)
                if length > 0:
                    for z in np.arange(0, length, 1):
                        # TODO: Review... no ho entenc (pq no afegir directament tot?) Ramon please check
                        self.PV_values.append(ray.PV_values[z])
            else:
                self.PV_energy.append(0.0)
                # TODO: Review... ray.wavelength always added to PV_wavelength
                self.PV_wavelength.append(ray.wavelength)
