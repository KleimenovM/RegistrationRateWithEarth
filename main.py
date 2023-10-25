# Evaluates average high-energy neutrino spectral attenuation
# for specific sources listed in "data/sources_table.csv"
import matplotlib.pyplot as plt
import numpy as np

from root_hist_draw import draw_root_hist
from single_theta_flux import calculate_single_theta_flux
from source import get_sources, Source
from telescope import Telescope, get_simple_telescope, get_Baikal, get_simple_telescope_from_txt
from transmission_function import TransmissionFunction


def get_simple_relative_flux(theta: np.ndarray, telescope: Telescope) -> np.ndarray:
    """
    Returns a relative flux value without attenuation in the Earth
    @param theta: np.ndarray - angular movement parametrization
    @param telescope: Telescope
    @return: np.ndarray - relative change in the flux due to the Earth revolution
    """
    lg_energy = telescope.lg_energy

    grid_x, grid_y = np.meshgrid(theta, lg_energy, indexing='ij')
    grid = np.array([grid_x, grid_y]).T

    registration_rate = telescope.ef_area(grid).T
    avg_reg_rate = np.mean(registration_rate, axis=0)

    return avg_reg_rate


def get_relative_flux(initial_flux: np.ndarray, theta: np.ndarray,
                      telescope: Telescope, tf: TransmissionFunction,
                      nuFate_method=1):
    """
    Returns e(mu) and tau relative fluxes with the initial flux given
    @param nuFate_method: 0 for no secondaries, 1 - e, mu with secondaries, 2 - tau with secondaries
    @param telescope: Telescope class object
    @param theta: trajectory parametrization with successive zenith angles
    @param initial_flux: flux on energy dependence before the Earth
    @param tf: Transmission Function class containing nuFate calculations
    @return: relative e(mu) flux, relative tau flux
    """
    m = theta.size
    n = telescope.lg_energy.size

    total_flux_matrix = np.zeros([m, n])

    mid_border, low_border = 1.0, -1.0

    for i, t_i in enumerate(theta):
        total_flux_matrix[i] \
            = calculate_single_theta_flux(initial_flux, t_i,
                                          telescope, tf, nuFate_method=nuFate_method,
                                          mid_border=mid_border, low_border=low_border)
    plt.show()

    total_flux = total_flux_matrix.mean(axis=0)

    return total_flux


def one_telescope_full_cycle(source: Source, tf: TransmissionFunction, telescope: Telescope,
                             simple: bool = False, angular_precision: int = 180) -> np.ndarray:
    """
    Performs a full calculations cycle for one source and one telescope
    @param source: Source - the neutrino source
    @param tf: TransmissionFunction - the Earth transmission function
    @param telescope: Telescope - the telescope for which calculation is performed
    @param simple: bool - if true no attenuation in the Earth is considered
    @param angular_precision: int - number of angles to describe the Earth's rotation
    @return: np.ndarray - registered spectrum (registration rate on energy)
    """
    zenith_angles = telescope.get_orbit_parametrization(source, angular_precision)[1]

    # print(f"{telescope.name}: {telescope.source_available_time(source)}")

    ref_energy = telescope.energy
    d_lg_e = telescope.lg_energy[1] - telescope.lg_energy[0]
    de = 10 ** (telescope.lg_energy + d_lg_e) - ref_energy

    ref_initial_flux = source.flux_on_energy_function(ref_energy)

    if simple:
        rel_flux = get_simple_relative_flux(zenith_angles, telescope)
    else:
        initial_flux = source.flux_on_energy_function(tf.energy)
        emu_at = get_relative_flux(initial_flux, zenith_angles, telescope, tf, nuFate_method=1)
        rel_flux = 1 * emu_at

    year_seconds = 3600 * 24 * 365
    multiplier = ref_initial_flux * de * year_seconds

    return rel_flux * multiplier


def main():
    # sources from file "source_table.csv" -- potential high-energy neutrino sources
    sources = get_sources("data/source_table.csv")

    # Baikal-GVD telescope 51°46′N 104°24'E
    # no dependence of effective area on the zenith angle
    baikal_simple = get_simple_telescope("BaikalGVD-trigger", latitude=[51, 46],
                                         filename="data/eff_area_single_cluster.root", histname="hnu_trigger",
                                         brd_angle=[30])

    # zenith-angle-dependent version
    baikal_complex = get_Baikal("data/eff_area_trigger", name_addition="")

    baikal = baikal_complex

    # KM3Net telescope 36°17'N 15°58'E
    # no dependence of effective area on the zenith angle (2023 data)
    km3net = get_simple_telescope_from_txt("KM3Net-trigger", latitude=[36, 16],
                                           filename="data/KM3Net-total.txt",
                                           brd_angle=[-12])

    # Earth transmission function calculated with nuFate
    tf = TransmissionFunction()

    source_numbers = [10]
    # source_numbers = [x for x in range(12)]

    value_s, value_c, value_km = 20*5, 20*5, 5

    baikal_r, km3net_r, baikal_s_r = [], [], []
    for sn in source_numbers:
        # source = Source(name="Test1", declination_angle=-np.pi/2, k0=1e-11, gamma=3, right_ascension_time=0.12)
        source = sources[sn]  # take one source from the list

        baikal_rel = one_telescope_full_cycle(source=source, tf=tf, telescope=baikal, simple=False)
        baikal_no_at_rel = one_telescope_full_cycle(source=source, tf=tf, telescope=baikal, simple=True)
        km3net_rel = one_telescope_full_cycle(source=source, tf=tf, telescope=km3net)

        baikal_r.append(baikal_rel)
        km3net_r.append(km3net_rel)
        baikal_s_r.append(baikal_no_at_rel)

        # i_s = np.round(np.sum(baikal_no_at_rel * value_s), 2)
        # i_c = np.round(np.sum(baikal_rel * value_c), 2)
        # i_km = np.round(np.sum(km3net_rel * value_km), 2)
        #
        # print(f'{source.name} & {i_s} & {i_c} & {np.round(i_c / i_s, 2)}'+r' \\')
        # print(f'{source.name} & {i_s} & {i_km} & {np.round(i_s / i_km, 2)}'+r' \\')

    draw_root_hist(sources=sources, source_numbers=source_numbers,
                   energy_c=baikal.energy, reg=baikal_r,
                   energy_s=km3net.energy, simple_reg=km3net_r,
                   value_c=20 * 5, value_s=5)
    return


if __name__ == '__main__':
    main()
