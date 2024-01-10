import numpy as np
import ROOT as rt

from transmission_function import TransmissionFunction
from single_theta_flux import convoluted_flux


def save_as_a_root_hist(filename: str, histname: str, x: np.ndarray, y: np.ndarray, z: np.ndarray, data: np.ndarray):
    hist = rt.TH3F(histname, histname,
                   x.size - 1, x,
                   y.size - 1, y,
                   z.size - 1, z)

    for i, x_i in enumerate(x):
        for j, y_j in enumerate(y):
            for k, z_k in enumerate(z):
                hist.Fill(x_i, y_j, z_k, data[i, j, k])

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
    angles = np.linspace(-np.pi/2, 0, m)  # angular distribution

    L = 20
    gammas = np.linspace(2.0, 4.0, L)  # gammas

    transmission_cof = np.zeros([L, m, n])

    for (i, gamma) in enumerate(gammas):
        gamma_flux = energy ** (-gamma)
        for (j, angle) in enumerate(angles):
            transmission_cof[i, j] = convoluted_flux(gamma_flux, angle, tf, nuFate_method=1)

    save_as_a_root_hist(filename="data/transmission_coefficient.root",
                        histname="tr_cof",
                        x=gammas, y=angles, z=tf.lg_energy,
                        data=transmission_cof)

    return


calculate_transmission_coefficients()
