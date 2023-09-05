# Evaluates average high-energy neutrino spectral attenuation
# for specific sources listed in "data/sources_table.csv"

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from transmission_function import TransmissionFunction
from source import Source, set_a_source


def read_source_data(number: int):
    """
    This function provides parameters of neutrino sources
    k0 (TeV-1 cm-2 s-1) [into standard * 1e4]
    gamma (no dim)
    e_cut (TeV)
    beta (no dim)
    """
    filename = "datasource_table.csv"
    data = pd.read_csv(filename, sep=",")

    d_n = data.T[number]

    return d_n.loc[['Source', 'delta', 'k0', 'gamma', 'e_cut', 'beta']]


def main():
    sources_list_file = "data/source_table.csv"
    tbl = read_sources_file(sources_list_file)

    tf = TransmissionFunction()

    return


if __name__ == '__main__':
    main()
