# Neutrino through Earth propagation
# Source class description

import numpy as np
from tools import deg_to_rad


class Source:
    """
    This class describes a typical stellar neutrino source
    """

    def __init__(self, name: str, declination_angle: float, k0: float, gamma: float, e_cut=None, beta=None):
        # position parameter
        self.name = name
        self.delta = declination_angle

        # spectrum parameters
        self.k0 = k0  # 1e-11 TeV-1 s-1 cm-2
        self.gamma = gamma  # spectrum attenuation coefficient

        # cut-off parameters
        self.e_cut = e_cut  # GeV
        self.beta = beta  # GeV

    def info(self):
        print(f"{self.name}, declination: {self.delta}")
        print(f"k0 = {self.k0} GeV-1 s-1 m-2, gamma = {self.gamma}")
        if self.e_cut and self.beta:
            print(f"e_cut: {self.e_cut} GeV, beta: {self.beta}")
        pass

    def flux_on_energy_function(self, energy):
        # exponential cutoff
        if self.e_cut != np.nan and self.beta:
            return self.k0 * (energy / 1000) ** (-self.gamma) * np.exp(-(energy / self.e_cut) ** self.beta)

        # just a polynomial spectrum
        return self.k0 * (energy / 1000) ** (-self.gamma)


def set_a_source(line) -> Source:
    return Source(name=line['Source'],
                  declination_angle=deg_to_rad([line['delta']]),  # deg to rad
                  k0=line['k0'] * 1e-11 * 1e4 * 1e-3,  # TeV-1 s-1 cm-2 to GeV-1 s-1 m-2
                  gamma=line['gamma'],
                  e_cut=line['e_cut'] * 1e3,  # TeV to GeV
                  beta=line['beta'])


if __name__ == '__main__':
    print("Not for direct use")
