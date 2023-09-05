# Neutrino through Earth propagation
# Telescope class description

from source import Source
from tools import *


class Telescope:
    """
    This class describes a typical neutrino telescope
    """
    def __init__(self, latitude: float):
        self.phi = latitude

    # reconstruction properties
    def get_orbit_parametrization(self, source: Source, psi: np.ndarray):
        vec = np.zeros([3, psi.size])
        vec[0], vec[1], vec[2] = sph_coord(r=1, theta=source.delta, phi=psi)
        rm = rot_matrix(self.phi)
        vec = np.dot(rm, vec)
        theta = np.arcsin(vec[2])
        return vec, theta


if __name__ == '__main__':
    print("Not for direct use")
