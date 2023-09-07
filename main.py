# Evaluates average high-energy neutrino spectral attenuation
# for specific sources listed in "data/sources_table.csv"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ROOT as rt
from scipy.interpolate import interp1d

from tools import deg_to_rad
from transmission_function import TransmissionFunction
from source import Source, set_a_source
from telescope import Telescope


def get_sources(filename: str) -> list[Source]:
    """
    This function provides parameters of neutrino sources
    k0 (TeV-1 cm-2 s-1) [into standard * 1e4]
    gamma (no dim)
    e_cut (TeV)
    beta (no dim)
    """
    data = pd.read_csv(filename, sep=',')

    sources = []
    for i in range(data.shape[0]):
        line_i = data.T[i].loc[['Source', 'delta', 'k0', 'gamma', 'e_cut', 'beta']]
        source_i: Source = set_a_source(line_i)
        sources.append(source_i)

    return sources


def get_telescope(name: str, declination: list, filename: str, histname: str = "hnu") -> Telescope:
    f = rt.TFile(filename, "read")
    hist = f.Get(histname)
    n = len(hist)

    # low_end, width, value
    data: np.ndarray = np.zeros([3, n])

    for i in range(n):
        data[0, i] = hist.GetBinLowEdge(i)  # low level
        data[1, i] = hist.GetBinWidth(i)  # bin width
        data[2, i] = hist.GetBinContent(i)  # bin average value

    return Telescope(name=name,
                     latitude=deg_to_rad(declination),
                     ef_area_table=data)


def get_relative_flux(theta: np.ndarray, initial_flux: np.ndarray, tf: TransmissionFunction):
    """

    @param theta:
    @param initial_flux:
    @param tf:
    @return:
    """
    m = theta.size
    n = tf.n

    # total flux matrix
    total_flux_matrix_emu = np.zeros([m, n])
    total_flux_matrix_tau = np.zeros([m, n])

    for i, t_i in enumerate(theta):
        print(np.rad2deg(t_i))
        if t_i > 0:  # if there's no attenuation
            total_flux_matrix_emu[i] = initial_flux
            total_flux_matrix_tau[i] = initial_flux
            continue
        total_flux_matrix_emu[i] = tf.convolution(t_i, initial_flux, nuFate_method=1)
        total_flux_matrix_tau[i] = tf.convolution(t_i, initial_flux, nuFate_method=2)

    total_flux_emu = total_flux_matrix_emu.mean(axis=0)
    total_flux_tau = total_flux_matrix_tau.mean(axis=0)

    return total_flux_emu / initial_flux, total_flux_tau / initial_flux


def attenuation_extrapolation(lg_energy: np.ndarray, spectra_ratio: np.ndarray, t00=1):
    """

    @param lg_energy:
    @param spectra_ratio:
    @param t00:
    @return:
    """
    # estimate derivatives
    t0, t1, t2 = lg_energy[0:3]
    f0, f1, f2 = np.log10(spectra_ratio[0:3])

    g0 = f0
    g1 = (f1 - f0) / (t1 - t0)
    g2 = (f2 - 2*f1 + f0) / (t1 - t0)**2

    # find a, b, c parameters
    matrix = np.array([[(t0 - t00)**2 * t0**2, (t0 - t00)**2 * t0, (t0 - t00)**2],
                       [2 * (t0 - t00) * (2 * t0 - t00) * t0, (t0 - t00) * (3*t0 - t00), 2*(t0 - t00)],
                       [3 * (2*t0 - t00)**2 - t00**2, 6 * t0 - t00, 2]])
    g_vector = np.array([g0, g1, g2])

    abc_vector = np.linalg.solve(matrix, g_vector)
    a, b, c = abc_vector

    # set extrapolation_function
    def extrapolation_function(e):
        t = np.log10(e)
        return 10**((t - t00)**2 * (a * t**2 + b * t + c))

    return extrapolation_function


def united_parts_interpolation(e1, e2, e3, f1, f2, f3):
    e = np.hstack([e1, e2, e3])
    f = np.hstack([f1, f2, f3])
    return interp1d(e, f)


def calculate_attenuation(source: Source, telescope: Telescope, tf: TransmissionFunction, angular_precision=180):
    """

    @param angular_precision:
    @param source:
    @param telescope:
    @param tf:
    @return:
    """
    print(source.info())

    # angles sample to describe the source's trajectory
    psi_sample = np.linspace(0, np.pi, angular_precision)
    vec, theta = telescope.get_orbit_parametrization(source, psi_sample)

    # standard energy sample do describe the spectrum
    e_sample = tf.energy
    de = tf.lg_energy[1] - tf.lg_energy[0]
    initial_flux = source.flux_on_energy_function(e_sample)

    # attenuated neutrino flux
    relative_flux_emu, relative_flux_tau = get_relative_flux(theta, initial_flux, tf)
    zero_border = 2
    e_mid = 10**(np.arange(zero_border, 3, de))
    mid_e_extrapolation_emu = attenuation_extrapolation(tf.lg_energy, relative_flux_emu, t00=zero_border)
    mid_e_extrapolation_tau = attenuation_extrapolation(tf.lg_energy, relative_flux_tau, t00=zero_border)

    # create a continuous flux on energy function
    e_low = 10**(np.arange(-1, zero_border, de))
    low_e_extrapolation_emu, low_e_extrapolation_tau = np.ones([e_low.size]), np.ones([e_low.size])

    full_transmission_spectrum_emu = united_parts_interpolation(e1=e_low, f1=low_e_extrapolation_emu,
                                                                e2=e_mid, f2=mid_e_extrapolation_emu(e_mid),
                                                                e3=e_sample, f3=relative_flux_emu)
    full_transmission_spectrum_tau = united_parts_interpolation(e1=e_low, f1=low_e_extrapolation_tau,
                                                                e2=e_mid, f2=mid_e_extrapolation_tau(e_mid),
                                                                e3=e_sample, f3=relative_flux_tau)

    full_e = 10**np.linspace(0, 10, 1000)

    plt.figure(figsize=(10, 6))
    plt.plot((1, 100), (1, 1), linestyle='dashed', color='black')
    plt.plot(full_e, full_transmission_spectrum_emu(full_e), label=r'$\nu_e$ or $\nu_\mu$')
    plt.plot(full_e, full_transmission_spectrum_tau(full_e), label=r'$\nu_\tau$')
    plt.scatter(e_sample, relative_flux_emu, 3)
    plt.scatter(e_sample, relative_flux_tau, 3)
    plt.scatter(e_mid, mid_e_extrapolation_emu(e_mid), 3)
    plt.scatter(e_mid, mid_e_extrapolation_tau(e_mid), 3)

    plt.fill_betweenx((-0.1, 1.2), 1e2, 1e3, color='red', alpha=.1)
    plt.fill_betweenx((-0.1, 1.2), 1e3, 1e10, color='blue', alpha=.1)
    plt.fill_betweenx((-0.1, 1.2), 1, 1e2, color='black', alpha=.2)

    plt.legend()
    plt.title(r'Затухание потока нейтрино при прохождении через Землю для $\Gamma = 1.36$')
    plt.grid(linestyle='dashed', color='gray')
    plt.xscale('log')
    plt.xlabel(r'$E_\nu,\ GeV$')
    plt.ylabel(r'$\Phi / \Phi_0$')

    plt.ylim(-0.1, 1.2)
    plt.xlim(1, 1e10)

    # plt.tight_layout()
    plt.show()

    return full_transmission_spectrum_emu, full_transmission_spectrum_tau


def main():
    # sources from file "source_table.csv" -- potential high-energy neutrino sources
    sources = get_sources("data/source_table.csv")

    # Baikal-GVD telescope 51◦46′N 104◦24'E
    telescope = get_telescope("Baikal", [51, 46], "data/eff_area.root")

    # Earth transmission function calculated with nuFate
    tf = TransmissionFunction()

    source_numbers = [5]

    initial, attenuated = [], []
    for sn in source_numbers:
        source = Source(name="Test1", declination_angle=-np.pi/2, k0=1, gamma=1.36)
        # source = sources[sn]
        initial_sn, attenuated_sn = calculate_attenuation(source, telescope, tf)
        initial.append(initial_sn)
        attenuated.append(attenuated_sn)

    return


if __name__ == '__main__':
    main()
