# Evaluates average high-energy neutrino spectral attenuation
# for specific sources listed in "data/sources_table.csv"
import matplotlib.pyplot as plt
import numpy as np

from root_hist_draw import draw_root_hist
from single_theta_flux import calculate_single_theta_flux
from source import get_sources
from telescope import Telescope, get_telescope
from transmission_function import TransmissionFunction


def get_relative_flux(initial_flux: np.ndarray, theta: np.ndarray,
                      telescope: Telescope, tf: TransmissionFunction):
    """
    Returns e(mu) and tau relative fluxes with the initial flux given
    @param telescope:
    @param theta: trajectory parametrization with successive zenith angles
    @param initial_flux: flux on energy dependence before the Earth
    @param tf: Transmission Function class containing nuFate calculations
    @return: relative e(mu) flux, relative tau flux
    """
    m = theta.size
    n = telescope.lg_energy.size

    # total flux matrix
    total_flux_matrix_emu = np.zeros([m, n])
    total_flux_matrix_tau = np.zeros([m, n])

    mid_border, low_border = 1.0, -1.0

    for i, t_i in enumerate(theta):
        total_flux_matrix_emu[i] \
            = calculate_single_theta_flux(initial_flux, t_i,
                                          telescope, tf, nuFate_method=1,
                                          mid_border=mid_border, low_border=low_border)
        total_flux_matrix_tau[i] \
            = calculate_single_theta_flux(initial_flux, t_i,
                                          telescope, tf, nuFate_method=2,
                                          mid_border=mid_border, low_border=low_border)

    plt.show()

    total_flux_emu = total_flux_matrix_emu.mean(axis=0)
    total_flux_tau = total_flux_matrix_tau.mean(axis=0)

    return total_flux_emu, total_flux_tau


def get_simple_relative_flux(theta: np.ndarray, telescope: Telescope):
    lg_energy = telescope.lg_energy

    grid_x, grid_y = np.meshgrid(theta, lg_energy, indexing='ij')
    grid = np.array([grid_x, grid_y]).T

    registration_rate = telescope.ef_area(grid).T
    avg_reg_rate = np.mean(registration_rate, axis=0)

    return avg_reg_rate


def main():
    # sources from file "source_table.csv" -- potential high-energy neutrino sources
    sources = get_sources("data/source_table.csv")

    # Baikal-GVD telescope 51◦46′N 104◦24'E
    telescope = get_telescope("Baikal", [51, 46], "data/eff_area.root")

    # Earth transmission function calculated with nuFate
    tf = TransmissionFunction()

    angular_precision = 180
    source_numbers = [10]

    ref_energy = telescope.energy
    d_lg_e = telescope.lg_energy[1] - telescope.lg_energy[0]
    de = 10 ** (telescope.lg_energy + d_lg_e) - ref_energy

    initial, simply_registered, registered = [], [], []
    for sn in source_numbers:
        # source = Source(name="Test1", declination_angle=-np.pi/2, k0=0.1, gamma=2)
        source = sources[sn]  # take one source from the list
        print(source.info())

        zenith_angles = telescope.get_orbit_parametrization(source, angular_precision)[1]
        initial_flux = source.flux_on_energy_function(tf.energy)

        ref_initial_flux = source.flux_on_energy_function(ref_energy)
        emu_at, tau_at = get_relative_flux(initial_flux, zenith_angles, telescope, tf)
        rel_flux_r = get_simple_relative_flux(zenith_angles, telescope)

        initial.append(ref_initial_flux)
        simply_registered.append(rel_flux_r)
        registered.append(2/3 * emu_at + 1/3 * tau_at)

    year_seconds = 3600 * 24 * 365

    print(telescope.lg_energy)

    for i in range(len(source_numbers)):
        draw_root_hist(sources, source_numbers, telescope.get_energy_bins(),
                       simply_registered * initial[i] * de * year_seconds, registered * initial[i] * de * year_seconds)
    return


if __name__ == '__main__':
    main()
