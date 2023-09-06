# Evaluates average high-energy neutrino spectral attenuation
# for specific sources listed in "data/sources_table.csv"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ROOT as rt

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
    hist = rt.TFile(filename, "read").Get(histname)
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
    d_theta = 1 / angular_precision

    # energy sample do describe the spectrum
    n = 200
    e_sample = 10**(np.linspace(3, 10, n))
    initial_flux = source.flux_on_energy_function(e_sample)

    plt.plot(e_sample, initial_flux)
    plt.xscale('log')
    plt.show()

    # start averaging
    total_flux_matrix_emu = np.zeros([angular_precision, n])
    total_flux_matrix_tau = np.zeros([angular_precision, n])

    for i, t_i in enumerate(theta):
        zenith_i = np.pi/2 - t_i
        if zenith_i < np.pi / 2:  # if there's no attenuation
            total_flux_matrix_emu[i] = initial_flux
            total_flux_matrix_tau[i] = initial_flux
            continue
        total_flux_matrix_emu[i] = tf.convolution(zenith_i, initial_flux, nuFate_method=1)
        total_flux_matrix_tau[i] = tf.convolution(zenith_i, initial_flux, nuFate_method=2)

    total_flux_emu = total_flux_matrix_emu.mean(axis=0)
    total_flux_tau = total_flux_matrix_tau.mean(axis=0)
    print(total_flux_emu.shape)

    return 0, 0


def main():
    # sources from file "source_table.csv" -- potential high-energy neutrino sources
    sources = get_sources("data/source_table.csv")

    # Baikal-GVD telescope 51◦46′N 104◦24'E
    telescope = get_telescope("Baikal", [51, 46], "data/eff_area.root")

    # Earth transmission function calculated with nuFate
    tf = TransmissionFunction()

    source_numbers = [1]

    initial, attenuated = [], []
    for sn in source_numbers:
        initial_sn, attenuated_sn = calculate_attenuation(sources[sn], telescope, tf)
        initial.append(initial_sn)
        attenuated.append(attenuated_sn)

    return


if __name__ == '__main__':
    main()
