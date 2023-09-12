import matplotlib.pyplot as plt
import numpy as np

from nuFATE.cascade_secs import get_eigs as get_eigs2
from tools import smart_division
from transmission_calculation import get_att_value_secs
from transmission_function import TransmissionFunction


def check_transmission_function():
    angles = np.linspace(-np.pi/2, 0, 180)
    energy = 10**np.linspace(3, 10, 200)

    gamma = 2
    flavor = 2
    spectrum = energy ** (-gamma)

    tf = TransmissionFunction()

    if_delta_marker = False

    w2, v2, ci2, energy_nodes2, phi_02 = get_eigs2(flavor, gamma,
                                                   "nuFATE/NuFATECrossSections.h5",
                                                   pure_spectrum=False)
    plt.figure(figsize=(10, 6))

    for j in range(0, len(angles), 5):
        a_j = angles[j]
        z_j = np.pi / 2 - angles[j]
        att_with_regeneration = get_att_value_secs(w2, v2, ci2,
                                                   energy_nodes2, z_j,
                                                   energy, phi_02, if_tau=False, absolute=if_delta_marker)
        att_with_regeneration_tau = get_att_value_secs(w2, v2, ci2,
                                                       energy_nodes2, z_j,
                                                       energy, phi_02, if_tau=True, absolute=if_delta_marker)

        att_tabular = smart_division(tf.convolution(a_j, spectrum, 1), spectrum)
        att_tabular_tau = smart_division(tf.convolution(a_j, spectrum, 2), spectrum)


        plt.subplot(2, 1, 1)

        plt.plot(energy, att_with_regeneration, color='royalblue', alpha=0.5)
        plt.plot(energy, att_with_regeneration_tau, color='purple', alpha=0.5)

        plt.plot(energy, att_tabular, color='red', alpha=.5)
        plt.plot(energy, att_tabular_tau, color='orange', alpha=.5)

        plt.subplot(2, 1, 2)
        plt.plot(att_tabular / att_with_regeneration, color='royalblue', alpha=.5)
        plt.plot(att_tabular_tau / att_with_regeneration_tau, color='red', alpha=.5)

    plt.subplot(2, 1, 1)
    plt.xscale('log')
    plt.subplot(2, 1, 2)
    plt.xscale('log')
    plt.show()

    return


if __name__ == '__main__':
    check_transmission_function()
