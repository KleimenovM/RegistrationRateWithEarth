# Earth Transmission Function (ETF) calculation
# Allowed energies from 10^3 to 10^10 GeV
# This file creates a data.npy file which contains points of the ETF

# Attention! It takes several minutes to execute this file!

import numpy as np
import ROOT as rt

from nuFATE.cascade import get_eigs as get_eigs1
from nuFATE.cascade_secs import get_eigs as get_eigs2
import nuFATE.earth as earth

from generate_energy_distribution import set_delta_function
from tools import extrapolating_spline

# Avogadro's number
N_A = 6.0221415e23


def get_att_value(w, v, ci, energy_nodes, zenith, E, phi_in, absolute=False):
    """
    This function calculates attenuation value with the use of NuFATE methods
    @param w: transmission equation eigenvalues - 1
    @param v: transmission equation eigenvalues - 2
    @param ci: transmission equation solution
    @param energy_nodes: energy bins log_10(GeV)
    @param zenith: zenith angle (rad, > pi/2)
    @param E: energy to interpolate
    @param phi_in: inner flux
    @param absolute: marker to fix whether the resulting flux is relative or absolute
    @return: attenuated flux (relative or absolute) at energy E
    """
    t = earth.get_t_earth(zenith) * N_A  # g/ cm^2
    if absolute:
        phi_sol = np.dot(v, (ci * np.exp(w * t))) * energy_nodes ** (-2)  # this is the attenuated flux
    else:
        phi_sol = np.dot(v, (ci * np.exp(w * t))) / phi_in  # this is phi/phi_initial, i.e. the relative attenuation
    return extrapolating_spline(E, energy_nodes, phi_sol, if_delta=absolute)


def get_att_value_secs(w, v, ci, energy_nodes, zenith, E, phi_in, if_tau=False, absolute=False):
    """
    This function calculates attenuation value with the use of NuFATE methods
    @param w: transmission equation eigenvalues - 1
    @param v: transmission equation eigenvalues - 2
    @param ci: transmission equation solution
    @param energy_nodes: energy bins log_10(GeV)
    @param zenith: zenith angle (rad, > pi/2)
    @param E: energy to interpolate
    @param phi_in: inner flux
    @param if_tau: marker to return the flux for the given particle or for tau-particle
    @param absolute: marker to fix whether the resulting flux is relative or absolute
    @return: attenuated flux (relative or absolute) at energy E
    """
    t = earth.get_t_earth(zenith) * N_A  # g / cm^2
    if absolute:
        e_nd = np.hstack([energy_nodes, energy_nodes])
        phi_sol = np.dot(v, (ci * np.exp(w * t))) * e_nd ** (-2)  # this is the attenuated flux
    else:
        phi_sol = np.dot(v, (ci * np.exp(w * t))) / phi_in  # this is phi/phi_0, i.e. the relative attenuation

    if if_tau:
        phi_sol1 = phi_sol[200:400]  # the tau bit.
        return extrapolating_spline(E, energy_nodes, phi_sol1, if_delta=absolute)

    phi_sol1 = phi_sol[0:200]  # the non-tau bit.
    return extrapolating_spline(E, energy_nodes, phi_sol1, if_delta=absolute)


def save_npy(filename: str, data: list[np.ndarray]):
    """
    Saves the attenuation matrix in a numpy file
    @param filename: "folder/filename.npy"
    @param data: list of 3D numpy arrays
    @return:
    """
    np.save("data/data_mod.npy", data)
    return


def save_root(filename, data, names):
    """
    Saves the attenuation matrix in a root file as a TH3F
    @param names: titles for histograms
    @param filename: "folder/filename.root"
    @param data: 3D numpy matrix
    @return:
    """

    file = rt.TFile.Open(filename, "RECREATE")
    for k in range(len(data)):
        file.WriteObject(data[k], names[k])

    file.Close()

    return


def calculate_transmission(if_npy: bool, if_root: bool):
    # Choose the flavor & index you want
    flavor = 2  # 1,2,3 = e, mu, tau; negative sign for antiparticles
    if_delta_marker = True

    # define possible angles
    z_angle_min, z_angle_max, m = np.pi/2, np.pi, 180
    z_angles = np.linspace(z_angle_min, z_angle_max, m)

    # define applicable energies
    lg_e_min, lg_e_max, n = 3, 10, 200
    lg_e = np.linspace(lg_e_min, lg_e_max, n)
    energy = 10**np.linspace(3, 10, n)

    root_hists, att_matrices = [], []  # attenuation matrices

    hist_names = ["No_regeneration", "With_regeneration", "tau_with_regeneration"]

    for k in range(3):
        # npy matrices
        att_matrices.append(np.zeros([m, n, n]))

        # root histograms
        hist_k = rt.TH3F(hist_names[k], hist_names[k],
                         m-1, z_angles,
                         n-1, lg_e,
                         n-1, lg_e)

        hist_k.GetXaxis().SetTitle("zenith angle, rad")
        hist_k.GetYaxis().SetTitle("lg(E_in / GeV)")
        hist_k.GetZaxis().SetTitle("lg(E_out / GeV)")

        root_hists.append(hist_k)

    for i in range(n):  # split by energy
        lg_e_i = lg_e[i]
        e_i = 10**lg_e_i
        # create a delta_function
        gamma = set_delta_function(e_i)

        # calculate eigenvalues
        # for simple cascades
        w1, v1, ci1, energy_nodes1, phi_01 = get_eigs1(flavor, gamma,
                                                       "nuFATE/NuFATECrossSections.h5",
                                                       pure_spectrum=True)
        # for cascades with secondaries
        w2, v2, ci2, energy_nodes2, phi_02 = get_eigs2(flavor, gamma,
                                                       "nuFATE/NuFATECrossSections.h5",
                                                       pure_spectrum=True)
        for j in range(m):  # split by angle
            if j % 30 == 0:
                print(f"{i}-{j}")
            z_j = z_angles[j]

            # NO SECONDARY PARTICLES
            att_no_regeneration = get_att_value(w1, v1, ci1,
                                                energy_nodes1, z_j,
                                                energy, phi_01, absolute=if_delta_marker)

            # INCLUDE SECONDARIES
            att_with_regeneration = get_att_value_secs(w2, v2, ci2,
                                                       energy_nodes2, z_j,
                                                       energy, phi_02, if_tau=False, absolute=if_delta_marker)
            att_with_regeneration_tau = get_att_value_secs(w2, v2, ci2,
                                                           energy_nodes2, z_j,
                                                           energy, phi_02, if_tau=True, absolute=if_delta_marker)

            # FILL THE MATRIX
            att_matrices[0][j, i] = att_no_regeneration
            att_matrices[1][j, i] = att_with_regeneration
            att_matrices[2][j, i] = att_with_regeneration_tau

            # FILL ROOT HISTOGRAMS
            for k in range(n):  # split by energy
                root_hists[0].Fill(z_j, lg_e[i], lg_e[k], att_no_regeneration[k])
                root_hists[1].Fill(z_j, lg_e[i], lg_e[k], att_with_regeneration[k])
                root_hists[2].Fill(z_j, lg_e[i], lg_e[k], att_with_regeneration_tau[k])

    if if_npy:
        # save as a .npy file
        save_npy("data/data_mod.npy", att_matrices)

    if if_root:
        # save as a .root file
        save_root("data/data_mod.root", root_hists, hist_names)

    return


if __name__ == '__main__':
    # Pay attention! It takes several minutes to execute the file!
    # Do not run it pointlessly
    calculate_transmission(if_root=True, if_npy=False)
