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
    # angles sample to describe the source's trajectory
    phi_sample = np.linspace(0, np.pi, angular_precision)

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
