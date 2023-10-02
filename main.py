# Evaluates average high-energy neutrino spectral attenuation
# for specific sources listed in "data/sources_table.csv"
import matplotlib.pyplot as plt
import numpy as np

from root_hist_draw import draw_root_hist
from single_theta_flux import calculate_single_theta_flux
from source import get_sources, Source
from telescope import Telescope, get_simple_telescope, get_complex_telescope
from transmission_function import TransmissionFunction


def get_Baikal(latitude=(51, 46)) -> Telescope:
    filenames = []
    angles = []
    p = 0.0
    filenames.append(f"data/eff_area/eff_area_{0.0}_{0.1}.root")
    angles.append(0.0)
    while p < 0.99:
        p_0_wr = np.round(p, 1)
        p_1_wr = np.round(p + .1, 1)
        filenames.append(f"data/eff_area/eff_area_{p_0_wr}_{p_1_wr}.root")
        angle_p = -np.mean((np.arcsin(p), np.arcsin(p + 0.1)))
        angles.append(angle_p)
        p += .1
    filenames.append(filenames[-1])
    angles.append(-np.pi/2)

    return get_complex_telescope(name='BaikalGVD',
                                 latitude=latitude,
                                 filenames=filenames,
                                 angles=angles)


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
    telescope1 = get_simple_telescope("BaikalGVD", latitude=[51, 46], filename="data/eff_area.root")
    telescope2 = get_Baikal()

    # Earth transmission function calculated with nuFate
    tf = TransmissionFunction()

    angular_precision = 180

    ref_energy = telescope1.energy
    d_lg_e = telescope1.lg_energy[1] - telescope1.lg_energy[0]
    de = 10 ** (telescope1.lg_energy + d_lg_e) - ref_energy

    year_seconds = 3600 * 24 * 365

    source_numbers = [11, 10]

    initial, double_simple_r, simply_registered, registered = [], [], [], []
    for sn in source_numbers:
        # source = Source(name="Test1", declination_angle=-np.pi/2, k0=0.1, gamma=3)
        source = sources[sn]  # take one source from the list
        print(source.info())
        print(telescope2.source_available_time(source))

        zenith_angles = telescope2.get_orbit_parametrization(source, angular_precision)[1]
        initial_flux = source.flux_on_energy_function(tf.energy)

        ref_initial_flux = source.flux_on_energy_function(ref_energy)
        emu_at1, tau_at1 = get_relative_flux(initial_flux, zenith_angles, telescope2, tf)
        rel_flux_r1 = get_simple_relative_flux(zenith_angles, telescope1)
        rel_flux_r2 = get_simple_relative_flux(zenith_angles, telescope2)

        initial.append(ref_initial_flux)
        multiplier = ref_initial_flux * de * year_seconds
        simply_registered.append(rel_flux_r2 * multiplier)
        double_simple_r.append(rel_flux_r1 * multiplier)
        registered.append((2/3 * emu_at1 + 1/3 * tau_at1) * multiplier)

    draw_root_hist(sources, source_numbers, telescope1.energy,
                   double_simple_r, simply_registered, registered)
    return


if __name__ == '__main__':
    main()
