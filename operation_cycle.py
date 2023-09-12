import numpy as np
from scipy.interpolate import interp1d

from telescope import Telescope
from transmission_function import TransmissionFunction


def attenuation_extrapolation(lg_energy: np.ndarray, spectra_ratio: np.ndarray, t00: float):
    """
    Extrapolates attenuated spectrum to the region where nuFate cannot perform calculations
    @param lg_energy: lg of neutrino energy
    @param spectra_ratio: relative attenuated spectrum
    @param t00: low-energy zero attenuation border
    @return: extrapolation function: e -> relative flux at (1e^{t00}, 1e3)
    """
    # estimate derivatives
    t0, t1, t2 = lg_energy[0:3]
    f0, f1, f2 = np.log10(spectra_ratio[0:3])

    g0 = f0
    g1 = (f1 - f0) / (t1 - t0)
    g2 = (f2 - 2*f1 + f0) / (t1 - t0)**2

    # find a, b, c parameters
    matrix = np.array([[(t0 - t00)**2 * t0**2, (t0 - t00)**2 * t0, (t0 - t00)**2],
                       [2 * (t0 - t00) * (2 * t0 - t00) * t0, (t0 - t00) * (3*t0 - t00), 2*(t0 - t00)],
                       [3 * (2*t0 - t00)**2 - t00**2, 6 * t0 - t00, 2]])
    g_vector = np.array([g0, g1, g2])

    abc_vector = np.linalg.solve(matrix, g_vector)
    a, b, c = abc_vector

    # set extrapolation_function
    def extrapolation_function(e):
        t = np.log10(e)
        return 10**((t - t00)**2 * (a * t**2 + b * t + c))

    return extrapolation_function


def united_parts_interpolation(e1, e2, e3, f1, f2, f3):
    """
    Unites three separate relative spectrum parts into one wide spectrum
    @param e1: low energy (zero attenuation)
    @param e2: middle energy (extrapolation region)
    @param e3: high energy (nuFate calculations)
    @param f1: identical unit value (no attenuation)
    @param f2: extrapolated flux
    @param f3: nuFate calculations
    @return: extrapolation function: e (in 1e-1, 1e10) -> relative flux
    """
    e = np.hstack([e1, e2, e3])
    f = np.hstack([f1, f2, f3])
    return interp1d(e, f)


def convoluted_flux(initial_flux: np.ndarray, zenith: float,
                    tf: TransmissionFunction, nuFate_method: int):
    """
    Step 1. Convolution -> Getting relative flux after transmission of the Earth
    @param zenith:
    @param initial_flux:
    @param tf:
    @param nuFate_method:
    @return:
    """
    convoluted = tf.convolution(zenith, initial_flux, nuFate_method=nuFate_method)
    return convoluted / initial_flux


def interpolated_flux(relative_flux: np.ndarray, tf: TransmissionFunction,
                      mid_border: float, low_border: float):
    """
    Step 2. Low energy extrapolation
    @param relative_flux:
    @param tf:
    @param mid_border:
    @param low_border:
    @return:
    """
    # Step 2. Extrapolation
    de = tf.lg_energy[1] - tf.lg_energy[0]
    e_mid = 10**(np.arange(mid_border, tf.energy[0], de))
    mid_e_extrapolation = attenuation_extrapolation(tf.lg_energy, relative_flux, t00=mid_border)

    e_low = 10**(np.arange(low_border, mid_border, de))
    low_e_extrapolation = np.ones([e_low.size])

    total_relative_spectrum = united_parts_interpolation(e1=e_low, f1=low_e_extrapolation,
                                                         e2=e_mid, f2=mid_e_extrapolation(e_mid),
                                                         e3=tf.energy, f3=relative_flux)

    return total_relative_spectrum


def eff_area_step(total_relative_spectrum_function, zenith: float, telescope: Telescope):
    """
    Step 3. Effective area consideration
    @param total_relative_spectrum_function:
    @param zenith:
    @param telescope:
    @return:
    """
    # 3.1. Data generation from interpolated values
    energy = telescope.energy
    relative_spectrum = total_relative_spectrum_function(energy)

    # 3.2. Multiplication with effective area
    grid_x, grid_y = np.meshgrid(zenith, energy, indexing='ij')
    grid = np.array([grid_x, grid_y]).T

    effective_area = telescope.ef_area(grid).T[0]

    return relative_spectrum * effective_area




def calculate_single_theta_flux(initial_flux: np.ndarray, zenith: float,
                                telescope: Telescope, tf: TransmissionFunction,
                                nuFate_method: int, mid_border: float = 1., low_border: float = -1.):
    conv_f = convoluted_flux(initial_flux, zenith, tf, nuFate_method)
    interp_f = interpolated_flux(conv_f, tf, mid_border, low_border)
    final_spectrum = eff_area_step(interp_f, zenith, telescope)
    return final_spectrum


if __name__ == '__main__':
    print("Not for direct use!")
