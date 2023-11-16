# Neutrino through Earth propagation
# Source class description

import numpy as np
import pandas as pd
from tools import deg_to_rad, galactic2equatorial, equatorial2galactic


class ExtendedSource:
    """
    This class describes a typical extended stellar neutrino source
    """

    def __init__(self, name: str, center_longitude, center_latitude,
                 k0, gamma, e_cut=None, beta=None, longitude_loc=None, latitude_loc=None):
        self.name = name

        # positioning parameters
        self.ll = center_longitude  # longitude (l), rad
        self.b = center_latitude  # latitude (b), rad

        # localization area
        if longitude_loc:
            self.l_loc = longitude_loc
        else:
            self.l_loc = 2 * np.pi
        if latitude_loc:
            self.b_loc = latitude_loc
        else:
            self.b_loc = 2 * np.pi

        # spectrum parameters
        self.k0 = k0  # 1e-11 TeV-1 s-1 cm-2, a function of l and b
        self.gamma = gamma  # spectrum attenuation coefficient, a function of l and b

        # cut-off parameters
        self.e_cut = e_cut  # GeV, a function of l and b
        self.beta = beta  # GeV, a function of l and b

    def info(self):
        print(f"{self.name}")
        print(f"pos = {galactic2equatorial(self.ll, self.b)}")
        pass

    def flux_on_energy_function(self, energy, ll, b):
        # exponential cutoff
        cut_off = 1  # a simple polynomial spectrum
        if self.e_cut != np.nan and self.beta:
            cut_off = np.exp(-(energy / self.e_cut(ll, b)) ** self.beta(ll, b))  # a high-energy cutoff

        return self.k0 * (energy / 1000) ** (-self.gamma) * cut_off


if __name__ == '__main__':
    print("Not for direct use")
