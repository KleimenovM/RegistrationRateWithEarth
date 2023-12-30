import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from root_hist_draw import draw_root_ext
from single_theta_flux import calculate_single_theta_flux
from source_extended import ExtendedSource, galactic_center
from telescope import Telescope, get_Baikal
from transmission_function import TransmissionFunction


def ext_source_cycle(source: ExtendedSource, tf: TransmissionFunction, telescope: Telescope,
                     angular_precision: int = 180):
    ref_energy = telescope.energy
    d_lg_e = telescope.lg_energy[1] - telescope.lg_energy[0]
    de = 10 ** (telescope.lg_energy + d_lg_e) - ref_energy

    ref_initial_flux_partition = np.zeros([source.num, telescope.energy.size])
    total_flux_partition = np.zeros([source.num, angular_precision, ref_energy.size])

    for i in range(source.num):
        if i % 10 == 0:
            print(i)
        ll_i, b_i = source.galactic_l_grid[i], source.galactic_b_grid[i]
        ref_initial_flux_partition[i, :] = source.flux_on_energy_function(ref_energy, ll=ll_i, b=b_i)
        initial_flux_i = source.flux_on_energy_function(tf.energy, ll=ll_i, b=b_i)

        dec_i, alpha_i = source.equatorial_dec_grid[i], source.equatorial_ra_grid[i]

        for j in range(angular_precision):
            zenith_ij = telescope.equatorial2telescope(dec_i, alpha_i + i / angular_precision * 2 * np.pi)

            total_flux_partition[i, j, :] = calculate_single_theta_flux(initial_flux=initial_flux_i,
                                                                        zenith=zenith_ij,
                                                                        telescope=telescope,
                                                                        tf=tf,
                                                                        nuFate_method=1)

    total_relative_flux_emu = np.mean(total_flux_partition, axis=1)

    total_rel = 1 * total_relative_flux_emu

    year_seconds = 3600 * 24 * 365
    multiplier = ref_initial_flux_partition * de * year_seconds

    return np.sum(total_rel * multiplier, axis=0)


def main():
    baikal_trig = get_Baikal("data/eff_area_5", name_addition="", histname="hnu_trigger")
    baikal_reco = get_Baikal("data/eff_area_trigger", name_addition="", histname="hnu_reco")

    # Earth transmission function calculated with nuFate
    tf = TransmissionFunction()

    # setting a source

    # sampling = [25, 50, 100, 200, 400, 800]
    # result1, result2 = [], []
    # for s in sampling:
    #     source = galactic_center(num=s)
    #     result1.append(ext_source_cycle(source=source, tf=tf, telescope=baikal_reco))
    #     result2.append(ext_source_cycle(source=source, tf=tf, telescope=baikal_trig))

    source = galactic_center(num=400)
    result1 = ext_source_cycle(source=source, tf=tf, telescope=baikal_reco)
    result2 = ext_source_cycle(source=source, tf=tf, telescope=baikal_trig)

    draw_root_ext(energy_s=baikal_trig.energy, energy_c=baikal_reco.energy,
                  simple_reg=result2, reg=result1, value_c=20 * 5, value_s=4 * 5, caption_pos='right')

    # plt.scatter(sampling, np.sum(result1, axis=1))
    # plt.scatter(sampling, np.sum(result2, axis=1))
    # plt.xscale('log')
    # plt.show()

    return


if __name__ == '__main__':
    main()
