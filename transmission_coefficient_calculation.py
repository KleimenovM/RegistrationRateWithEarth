import numpy as np
import ROOT as rt

from transmission_function import TransmissionFunction
from single_theta_flux import convoluted_flux


def save_as_a_2d_root_hist(filename: str, histname: str, x: np.ndarray, y: np.ndarray, data: np.ndarray,
                           x_title: str = "zenith angle, rad", y_title: str = "lg(E/1 GeV)"):
    hist = rt.TH2F(histname, histname,
                   x.size - 1, x,
                   y.size - 1, y)

    for i, x_i in enumerate(x):
        for j, y_j in enumerate(y):
            hist.Fill(x_i, y_j, data[i, j])

    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)

    file = rt.TFile.Open(filename, "RECREATE")
    file.WriteObject(hist, histname)

    file.Close()

    return


def calculate_transmission_coefficients():
    """

    @return:
    """
    tf = TransmissionFunction()

    n = tf.n
    energy = tf.energy  # energy distribution (from 10^3 GeV to 10^10 GeV, 200 points)

    m = 180
    angles = np.linspace(np.pi/2, np.pi, m)  # angular distribution

    L = 10
    gammas = np.linspace(2.0, 4.0, L)  # gammas

    transmission_cof = np.zeros([L, m, n])

    for (i, gamma) in enumerate(gammas):
        gamma_flux = energy ** (-gamma)

        print(gamma)

        for (j, angle) in enumerate(angles):

            transmission_cof[i, j] = convoluted_flux(gamma_flux, np.pi/2 - angle, tf, nuFate_method=1)

        save_as_a_2d_root_hist(filename=f"data/transmission/tc_gamma{np.round(gamma, 1)}.root",
                               histname=f"tc",
                               x=angles, y=tf.lg_energy,
                               data=transmission_cof[i])

    return


calculate_transmission_coefficients()
