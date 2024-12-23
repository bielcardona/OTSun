"""Module otsun.experiments for creating experiments.

The module defines a class `Experiment` that deals with the setup and
run of experiments.
"""
from .movements import MultiTracking
from .source import LightSource, GeneralizedSunWindow, buie_distribution
from .scene import Scene
from .math import polar_to_cartesian, cdf_from_pdf_file
from FreeCAD import Document
from importlib.resources import files

import numpy as np
import logging

logger = logging.getLogger(__name__)


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
    document_to_show : App.Document
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

    def __init__(
            self,
            scene: Scene,
            light_source: LightSource,
            number_of_rays: int,
            document_to_show: Document | None = None
    ):
        self.scene = scene
        self.light_source = light_source
        self.document_to_show = document_to_show
        if document_to_show:
            self.light_source.emitting_region.add_to_document(document_to_show)
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

    def run(self) -> None:
        """Runs the experiment and plots the rays in the document specified (if any)"""
        for n in np.arange(0, self.number_of_rays, 1):
            logger.info(f"Emitting ray {n}")
            ray = self.light_source.emit_ray()
            ray.run()
            self.wavelengths.append(ray.wavelength)
            if self.document_to_show:
                ray.add_to_document(self.document_to_show)
            if ray.Th_absorbed:
                self.captured_energy_Th += ray.energy
                self.Th_energy.append(ray.energy)
                self.Th_wavelength.append(ray.wavelength)
                self.points_absorber_Th.append(
                    (ray.energy,
                     ray.points[-1].x, ray.points[-1].y, ray.points[-1].z,
                     ray.points[-2].x, ray.points[-2].y, ray.points[-2].z,
                     ray.last_normal.x, ray.last_normal.y, ray.last_normal.z,
                     ray.wavelength)
                )
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


class TotalAnalysis:
    def __init__(
            self,
            scene: Scene,
            number_of_rays: int,
            phi: float,
            theta: float,
            track_movements: bool = False,
            csr_value: float | None = None,
            aperture_th: float | None = None,
            aperture_pv: float | None = None,
    ):
        self.scene = scene
        self.number_of_rays = number_of_rays
        self.phi = phi
        self.theta = theta
        self.csr_value = csr_value
        self.track_movements = track_movements
        self.has_run = False
        self.aperture_th = aperture_th
        self.aperture_pv = aperture_pv

        self.main_direction = polar_to_cartesian(phi, theta) * (-1.0)
        if self.track_movements:
            tracking = MultiTracking(self.main_direction, self.scene)
            tracking.make_movements()
        emitting_region = GeneralizedSunWindow(self.scene, self.main_direction)
        self.aperture = emitting_region.aperture
        data_file_spectrum = files('otsun').joinpath('data', 'ASTMG173-direct.txt').__str__()
        light_spectrum = cdf_from_pdf_file(data_file_spectrum)
        if csr_value is None:
            direction_distribution = None
        else:
            direction_distribution = buie_distribution(csr_value)
        light_source = LightSource(scene, emitting_region, light_spectrum, 1.0,
                                   direction_distribution, None)
        self.experiment = Experiment(scene, light_source, number_of_rays, None)

    def run(self):
        self.experiment.run()
        self.has_run = True

    @property
    def results(self):
        if not self.has_run:
            raise ValueError("Experiment has not been run")
        if self.aperture_pv is not None and self.aperture_pv != 0:
            efficiency_pv = (
                    (self.experiment.captured_energy_PV / self.aperture_pv) /
                    (self.number_of_rays / self.aperture))
        else:
            efficiency_pv = 0
        if self.aperture_th is not None and self.aperture_th != 0:
            efficiency_th = (
                    (self.experiment.captured_energy_Th / self.aperture_th) /
                    (self.number_of_rays / self.aperture))
        else:
            efficiency_th = 0

        return {
            'captured_Th': self.experiment.captured_energy_Th,
            'captured_PV': self.experiment.captured_energy_PV,
            'aperture': self.experiment.light_source.emitting_region.aperture,
            'efficiency_Th': efficiency_th,
            'efficiency_PV': efficiency_pv,
        }
