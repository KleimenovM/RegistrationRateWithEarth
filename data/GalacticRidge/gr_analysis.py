import numpy as np
import matplotlib.pyplot as plt


def read_data(filename):
    lg_energy, e2_flux = np.loadtxt(filename, unpack=True, skiprows=2)
    energy = 10 ** lg_energy  # GeV
    flux = e2_flux / energy ** 2
    return energy, flux
