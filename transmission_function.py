import numpy as np
from scipy.interpolate import RegularGridInterpolator as interp3d
from scipy.integrate import trapezoid, simpson
import matplotlib.pyplot as plt


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
        theta_min, theta_max, self.m = 0, -np.pi/2, 180
        lg_e_min, lg_e_max, self.n = 3, 10, 200
        # angle and log_energy axis
        self.angles = np.linspace(theta_min, theta_max, self.m)
        self.lg_energy = np.linspace(lg_e_min, lg_e_max, self.n)
        self.energy = 10 ** self.lg_energy
        # regeneration functions
        self.no_regen_function = None
        self.with_regen_function = None
        self.tau_regen_function = None
        # interpolation
        self.interpolated_table()

    def interpolated_table(self):
        """
        Interpolates tabular data (calculated with nuFate)
        """
        xyz = (self.angles, self.lg_energy, self.lg_energy)

        no_regeneration = self.table_data[0]
        with_regeneration = self.table_data[1]
        tau_regeneration = self.table_data[2]

        self.no_regen_function = interp3d(xyz, no_regeneration, method="linear")
        self.with_regen_function = interp3d(xyz, with_regeneration, method="linear")
        self.tau_regen_function = interp3d(xyz, tau_regeneration, method="linear")
        pass

    def convolution(self, angle, input_spectrum, integration_method=trapezoid, nuFate_method=0):
        """
        Performs convolution of the input spectrum with the Earth's transmission function
        @param angle: zenith angle
        @param input_spectrum: numpy array (self.n) - flux on angle & energy dependence
        @param integration_method: function that takes (y, x, axis=0) and returns an integral
        @param nuFate_method: 0 -> cascades only, 1 -> secondaries included, 2-> tau neutrinos
        @return: numpy array (self.n) - attenuated flux on energy dependence
        """
        transmission_functions = (self.no_regen_function, self.with_regen_function, self.tau_regen_function)
        tf = transmission_functions[nuFate_method]

        # calculate final flux for each angle
        grid_x, grid_y, grid_z = np.meshgrid(angle, self.lg_energy, self.lg_energy, indexing='ij')
        grid = np.array([grid_x, grid_y, grid_z]).T

        angle_transmission = tf(grid).T[0]

        input_spectrum_matrix_a = np.repeat(input_spectrum, self.n).reshape(self.n, self.n).T
        product = np.dot(input_spectrum_matrix_a, angle_transmission)

        return integration_method(product, self.energy, axis=0) / (self.energy[-1] - self.energy[0])


if __name__ == '__main__':
    transmission_function = TransmissionFunction()
    e = transmission_function.lg_energy
    de = e[1] - e[0]

    ang = transmission_function.angles[10]
    energy_flux_matrix = (10 ** e) ** (-1.36)

    final = transmission_function.convolution(ang, energy_flux_matrix)

    plt.plot(10**e, final, label='final flux')
    plt.plot(10**e, energy_flux_matrix, label='initial flux')
    plt.plot(10**e, final / energy_flux_matrix, label='ratio')
    plt.xscale('log')
    # plt.yscale('log')
    plt.legend()
    plt.show()
