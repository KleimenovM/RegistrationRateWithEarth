import numpy as np
from scipy.interpolate import RegularGridInterpolator as interp3d
from scipy.integrate import trapezoid, simpson
import matplotlib.pyplot as plt

from tools import smart_division


class TransmissionFunction:
    """
    This class takes a matrix, calculated with nuFate, from data_mod.npy
    and interpolates it in order to simplify the usage of nuFate package

    To use the interpolation, there is the -convolution- function which
    results into attenuated flux given an initial flux and zenith angle
    """

    def __init__(self):
        # load data.npy file
        self.table_data = np.load("data/data_mod.npy")
        # energies and angles boundaries
        theta_min, theta_max, self.m = 0, -np.pi / 2, 180
        lg_e_min, lg_e_max, self.n = 3, 10, 200
        # angle and log_energy axis
        self.angles = np.linspace(theta_min, theta_max, self.m)
        self.lg_energy = np.linspace(lg_e_min, lg_e_max, self.n)
        self.energy = 10 ** self.lg_energy

    def angle_interpolated_matrix(self, angle: float, nuFate_method: int):
        if angle >= 0:
            return np.eye(self.n)  # no attenuation

        transmission_table = self.table_data[nuFate_method]

        j = self.angles[self.angles > angle].size - 1
        t = (angle - self.angles[j]) / (self.angles[j + 1] - self.angles[j])

        transmission_matrix = transmission_table[j] * (1 - t) + transmission_table[j + 1] * t
        return transmission_matrix

    def convolution(self, angle, input_spectrum, nuFate_method: int):
        """
        Performs convolution of the input spectrum with the Earth's transmission function
        @param angle: zenith angle
        @param input_spectrum: numpy array (self.n) - flux on angle & energy dependence
        @param nuFate_method: 0 -> cascades only, 1 -> secondaries included, 2-> tau neutrinos
        @return: numpy array (self.n) - attenuated flux on energy dependence
        """
        if angle >= 0:
            return input_spectrum

        angle_transmission = self.angle_interpolated_matrix(angle, nuFate_method)
        product = np.dot(input_spectrum, angle_transmission)

        return product


if __name__ == '__main__':
    print("Not for direct use")
