# Neutrino through Earth propagation
# Telescope class description
import numpy as np
from scipy.interpolate import interp1d

from source import Source
from tools import sph_coord, rot_matrix


class Telescope:
    """
    This class describes a typical neutrino telescope
    """
    def __init__(self, name: str, latitude: float, ef_area_table: np.ndarray):
        self.name = name
        self.phi = latitude
        self.ef_area_table = ef_area_table
        self.ef_area = None
        self.set_ef_area()

    # reconstruction properties

    def set_ef_area(self):
        e = self.ef_area_table[0] + self.ef_area_table[1] / 2  # middle of the bin
        value = self.ef_area_table[2]  # ef_area value
        self.ef_area = interp1d(e, value, kind='cubic')
        pass

    def get_orbit_parametrization(self, source: Source, psi: np.ndarray):
        vec = np.zeros([3, psi.size])
        vec[0], vec[1], vec[2] = sph_coord(r=1, theta=source.delta, phi=psi)
        rm = rot_matrix(self.phi)
        vec = np.dot(rm, vec)
        theta = np.arcsin(vec[2])
        return vec, theta


if __name__ == '__main__':
    print("Not for direct use")
