# Neutrino through Earth propagation
# Source class description

import numpy as np
import pandas as pd
from tools import galactic2equatorial, equatorial2galactic, get_a_const_func


class ExtendedSource:
    """
    This class describes a typical extended stellar neutrino source
    """

    def __init__(self, name: str, center_longitude, center_latitude,
                 k0, gamma, e_cut=None, beta=None, longitude_loc=None, latitude_loc=None,
                 num=100):
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

        self.num = num

        gal_grid, eq_grid = self.single_rectangular_sampling()
        self.galactic_l_grid, self.galactic_b_grid = gal_grid[0], gal_grid[1]
        self.equatorial_dec_grid, self.equatorial_ra_grid = eq_grid[0], eq_grid[1]

    def info(self):
        print(f"{self.name}")
        print(f"pos = {galactic2equatorial(self.ll, self.b)}")
        pass

    def flux_on_energy_function(self, energy, ll, b):
        # exponential cutoff
        cut_off = 1  # a simple polynomial spectrum
        if self.e_cut != np.nan and self.beta:
            cut_off = np.exp(-(energy / self.e_cut(ll, b)) ** self.beta(ll, b))  # a high-energy cutoff

        delta_sin = 2 * np.sin(self.b_loc/2)

        return (self.k0 * (energy / 1000) ** (-self.gamma) * cut_off) / self.num

    def single_rectangular_sampling(self):
        total_area = self.l_loc * self.b_loc  # rad^2
        single_cell_area = total_area / self.num  # rad^2
        cell_size = np.sqrt(single_cell_area)  # rad, angular size of a squared cell

        nx = int(self.l_loc // cell_size)
        ny = int(self.b_loc // cell_size)

        self.num = nx * ny

        l_split = self.ll + np.linspace(-self.l_loc / 2, self.l_loc / 2, nx)
        b_split = self.b + np.linspace(-self.b_loc / 2, self.b_loc / 2, ny)

        gal_ll_grid, gal_b_grid = np.zeros(self.num), np.zeros(self.num)
        eq_dec_grid, eq_ra_grid = np.zeros(self.num), np.zeros(self.num)

        for i in range(nx):
            for j in range(ny):
                gal_ll_grid[i * ny + j], gal_b_grid[i * ny + j] = l_split[i], b_split[j]
                eq_dec_grid[i * ny + j], eq_ra_grid[i * ny + j] = galactic2equatorial(
                    l_split[i], b_split[j], get_delta=True, get_alpha=True)

        return [gal_ll_grid, gal_b_grid], [eq_dec_grid, eq_ra_grid]


def galactic_center(num=100):
    # Data from this article
    # https://arxiv.org/pdf/2307.01038.pdf
    ll, b = .0, .0  # rad
    l_loc, b_loc = np.deg2rad(60), np.deg2rad(4)  # [-30, 30] & [-2, 2]
    k0_ref = 1.66 * 1e-18  # GeV-1 cm-2 s-1 sr-1
    gamma_val = 2.53
    k0_val = k0_ref * 1e4 * (100**gamma_val)
    k0_f = get_a_const_func(k0_val)
    gamma_f = get_a_const_func(gamma_val)

    return ExtendedSource(name='Galactic Center',
                          center_longitude=ll, center_latitude=b,
                          k0=k0_f, gamma=gamma_f,
                          longitude_loc=l_loc, latitude_loc=b_loc,
                          num=num)


if __name__ == '__main__':
    print("Not for direct use")
