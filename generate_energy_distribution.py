import numpy as np
import matplotlib.pyplot as plt


Lg_min, Lg_max, N = 3, 10, 200
D_lg_e = (Lg_max - Lg_min) / N
Lg_e = np.linspace(3, 10, 200)


def set_delta_function(e):
    """
    Returns a delta-spectrum with peak at e
    @param e: delta function parameter
    @return: np.ndarray -- delta spectrum (_____|_______)
    """
    log10_e = np.log10(e)
    delta_position = np.all([Lg_e - D_lg_e / 2 <= log10_e, Lg_e + D_lg_e / 2 > log10_e], axis=0)
    e_bin = np.argmax(delta_position)
    n = Lg_e.size
    # if e_bin == n - 1:
    #     width = 10**(Lg_e[e_bin] + D_lg_e) - 10**Lg_e[e_bin]
    # else:
    #     width = 10**Lg_e[e_bin + 1] - 10**Lg_e[e_bin]
    width = 1
    return delta_position / width


def draw_a_figure(energy, source_flux):
    """
    Draws a given flux on log_10(energy)
    @param energy: x-axis (energy)
    @param source_flux: y-axis (flux)
    @return: nothing
    """
    plt.plot(energy, 10 ** Lg_e, source_flux)
    plt.xscale('log')
    plt.show()
    return


if __name__ == '__main__':
    # Draw a delta-function
    n = 20
    energy_i = 10**Lg_e[n]  # GeV
    flux = set_delta_function(energy_i)
    draw_a_figure(energy_i, flux)
