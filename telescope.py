# Neutrino through Earth propagation
# Telescope class description
import numpy as np
from scipy.interpolate import RegularGridInterpolator as interp2d

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
        self.lg_energy = None
        self.energy = None
        self.ef_area = None
        self.simple_effective_area()

    def get_orbit_parametrization(self, source: Source, angular_precision: int):
        psi = np.linspace(0, 2 * np.pi, angular_precision)
        vec = np.zeros([3, psi.size])
        vec[0], vec[1], vec[2] = sph_coord(r=1, theta=source.delta, phi=psi)
        rm = rot_matrix(self.phi)
        vec = np.dot(rm, vec)
        theta = np.arcsin(vec[2])
        return vec, theta

    def source_available_time(self, source):
        m = 1000
        vec, theta = self.get_orbit_parametrization(source, m)
        theta_good = theta < - np.pi / 6
        return np.sum(theta_good) / m

    def simple_effective_area(self, angle_precision: int = 180):
        self.lg_energy = self.ef_area_table[0]  # + self.ef_area_table[1] / 2  # middle of the bin
        self.energy = 10**self.lg_energy
        value = self.ef_area_table[2]  # ef_area value

        zenith_parametrization = np.linspace(-np.pi/2, np.pi/2, angle_precision)
        ef_area_parametrization = np.zeros([angle_precision, value.size])

        # simple method: if zenith angle < border angle, ef. area = 0
        border_angle = -np.pi/6

        for i, z in enumerate(zenith_parametrization):
            if z < border_angle:
                ef_area_parametrization[i] = value

        xy = (zenith_parametrization, self.lg_energy)

        self.ef_area = interp2d(xy, ef_area_parametrization, method='linear')
        pass


if __name__ == '__main__':
    print("Not for direct use")
