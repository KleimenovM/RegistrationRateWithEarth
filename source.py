# Neutrino through Earth propagation
# Source class description


from tools import *


class Source:
    """
    This class describes a typical stellar neutrino source
    """
    def __init__(self, declination_angle: float, k0: float, gamma: float, e_cut=None, beta=None):
        # position parameter
        self.delta = declination_angle

        # spectrum parameters
        self.k0 = k0  # 1e-11 TeV-1 s-1 cm-2
        self.gamma = gamma  # spectrum attenuation coefficient

        # cut-off parameters
        self.e_cut = e_cut  # GeV
        self.beta = beta  # GeV

    def flux_on_energy_function(self, energy):
        # exponential cutoff
        if self.e_cut and self.beta:
            return self.k0 * (energy / 1000)**(-self.gamma) * np.exp(-(energy / self.e_cut)**self.beta)

        # just a polynomial spectrum
        return self.k0 * (energy / 1000)**(-self.gamma)


def set_a_source(line: np.ndarray):
    return


if __name__ == '__main__':
    print("Not for direct use")
