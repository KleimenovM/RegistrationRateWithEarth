import numpy as np
import ROOT as rt
import matplotlib.pyplot as plt

from telescope import Telescope, get_Baikal


def check_ef_area():

    baikal: Telescope = get_Baikal("data/eff_area_trigger")

    lg_energy_range = baikal.lg_energy
    energy_range = baikal.energy

    a = -np.pi/3

    xv = np.array(np.meshgrid(a, lg_energy_range, indexing='ij')).T
    f_xv = baikal.ef_area(xv).T[0]

    plt.scatter(energy_range, f_xv)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

    return


if __name__ == '__main__':
    check_ef_area()
