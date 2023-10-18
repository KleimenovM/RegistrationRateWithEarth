import numpy as np
import matplotlib.pyplot as plt

from transmission_function import TransmissionFunction
from single_theta_flux import convoluted_flux, interpolated_flux


def extrapolation_check():

    tf = TransmissionFunction()

    angle = -np.pi / 3  # zenith -60 deg
    lg_e0 = np.linspace(3, 10, 200)
    e0 = 10**lg_e0
    gamma = 2

    initial_flux = e0 ** (-gamma)

    conv_f_emu = convoluted_flux(initial_flux, angle, tf, nuFate_method=1)
    conf_f_tau = convoluted_flux(initial_flux, angle, tf, nuFate_method=2)

    tr_flux_emu = interpolated_flux(conv_f_emu, tf, mid_border=2, low_border=-1)
    tr_flux_tau = interpolated_flux(conf_f_tau, tf, mid_border=2, low_border=-1)

    lg_e = np.linspace(0, 10, 200)
    e = 10**lg_e

    plt.figure(figsize=(8, 6))

    plt.plot(e, tr_flux_emu(e), linewidth=1.5, color="royalblue", label=r"$\nu_e\ or\ \nu_\mu$")
    plt.plot(e, tr_flux_tau(e), linewidth=1.5, color="#c20", label=r"$\nu_\tau$")

    plt.fill_betweenx((-0.2, 1.5), (1e0, 1e0), (10**1.5, 10**1.5), color='gray', alpha=.2, linewidth=0)
    plt.fill_betweenx((-0.2, 1.5), (10**1.5, 10**1.5), (1e3, 1e3), color='red', alpha=.2, linewidth=0)
    plt.fill_betweenx((-0.2, 1.5), (1e3, 1e3), (1e10, 1e10), color='blue', alpha=.2, linewidth=0)

    plt.xscale('log')

    plt.xlim(1, 1e10)
    plt.ylim(-0.1, 1.2)
    plt.xlabel(r"$E,\ GeV$", fontsize=14)
    plt.ylabel(r"$\phi/\phi_0$", fontsize=14)

    plt.tick_params(labelsize=14)

    plt.legend(fontsize=14)

    plt.grid(linestyle='dashed', color='gray')

    plt.tight_layout()
    plt.show()

    return


if __name__ == '__main__':
    extrapolation_check()
