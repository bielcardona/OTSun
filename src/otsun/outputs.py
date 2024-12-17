"""Module otsun.outputs

Helper functions to format data for output
"""
import os

import numpy as np

def spectrum_to_constant_step(
        file_in: str | bytes | os.PathLike,
        wavelength_step: float,
        wavelength_min: float,
        wavelength_max: float
) -> np.ndarray:
    data_array = np.loadtxt(file_in, usecols=(0, 1))
    wl_spectrum = data_array[:, 0]
    I_spectrum = data_array[:, 1]
    array_inter = [[x, np.interp(x, wl_spectrum, I_spectrum)] for x in
                   np.arange(wavelength_min, wavelength_max + wavelength_step / 2.0, wavelength_step)]
    return np.array(array_inter)


def make_histogram_from_experiment_results(
        results_wavelength : list[np.ndarray],
        results_energy: list[np.ndarray],
        step_wavelength: float,
        aperture_collector: float,
        aperture_source: float
) -> np.ndarray:
    y_energy = np.array(np.concatenate(results_energy))
    y_energy = (y_energy / aperture_collector) / (1.0 / aperture_source)
    x_wavelength = np.array(np.concatenate(results_wavelength))
    data_ = np.array([x_wavelength, y_energy])
    data_ = data_.T
    bins_ = np.arange(int(np.amin(x_wavelength)), np.amax(x_wavelength) + step_wavelength * 1.1, step_wavelength)
    table_ = np.histogram(data_[:, 0], bins=bins_, weights=data_[:, 1])
    norm = np.histogram(data_[:, 0], bins=bins_)
    u = np.divide(table_[0], norm[0])
    bins = np.arange(int(np.amin(x_wavelength)), np.amax(x_wavelength) + step_wavelength, step_wavelength)
    table_ = np.column_stack((bins, u))
    return table_


def twoD_array_to_constant_step(
        twoD_array: np.ndarray,
        step: float,
        wavelength_min: float,
        wavelength_max: float
) -> np.ndarray:
    wl_spectrum = twoD_array[:, 0]
    I_spectrum = twoD_array[:, 1]
    array_inter = [[x, np.interp(x, wl_spectrum, I_spectrum)] for x in
                   np.arange(wavelength_min, wavelength_max + step / 2.0, step)]
    return np.array(array_inter)


def spectral_response(
        optical_absorption_wavelength: np.ndarray,
        iqe: float | str | bytes | os.PathLike
) -> np.ndarray:
    q_e = 1.60217662E-19
    h = 6.62607E-34
    c = 299792458.0
    hc = h * c
    opt_wavelength = optical_absorption_wavelength
    if np.isreal(iqe):
        SR = [[opt[0], iqe * opt[0] * opt[1] * q_e * 1E-9 / hc, ] for opt in opt_wavelength]
    else:
        data_array = np.loadtxt(iqe, usecols=(0, 1))
        wl_ = data_array[:, 0]
        iqe_ = data_array[:, 1]
        SR = [[opt[0], np.interp(opt[0], wl_, iqe_) * opt[0] * opt[1] * q_e * 1E-9 / hc, ] for opt in opt_wavelength]
    return np.array(SR)


def photo_current(
        spectral_response: np.ndarray,
        source_spectrum: np.ndarray
) -> float:
    wavelengths = source_spectrum[:, 0]
    SR_by_spectrum = spectral_response[:, 1] * source_spectrum[:, 1]
    photo_current = np.trapz(SR_by_spectrum, x=wavelengths)
    return photo_current


# ---
# Helper functions for outputs in Total Analysis
# ---

def integral_from_data_file(file_in):
    source_spectrum = np.loadtxt(file_in, usecols=(0, 1))
    integral = np.trapz(source_spectrum[:, 1], x=source_spectrum[:, 0])
    return integral
