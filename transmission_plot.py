import matplotlib.pyplot as plt
import numpy as np
import ROOT as rt

from tools import smart_division
from transmission_function import TransmissionFunction


def plot_transmission_graph(gamma=2):
    tf = TransmissionFunction()

    angles = np.linspace(-np.pi / 2, 0, 180)
    e = tf.energy

    flux = e ** (-gamma)

    m, n = angles.size, e.size

    hist = rt.TH2F("Title", "Transmission matrix for #Gamma = 2", n - 1, e, m - 1, np.rad2deg(angles))

    for i, a in enumerate(angles):
        f_i = smart_division(tf.convolution(a, flux, 2), flux)
        for j, e_j in enumerate(e):
            hist.Fill(e_j, np.rad2deg(a), f_i[j])

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.1)
    canvas.SetBottomMargin(.1)
    canvas.SetRightMargin(.18)

    size = 0.04

    axes = hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()
    labels = ["E, GeV", "#theta, deg", "#phi/#phi_{0}"]

    for i, axis in enumerate(axes):
        axis.SetTitle(labels[i])
        axis.SetTitleSize(size)
        axis.SetLabelSize(size)
        axis.SetTitleOffset(1.2)

    # Y-AXIS LABELLING SHIFT
    y_axis = axes[1]
    for i in range(10):
        y_axis.ChangeLabel(i+1, -1, -1, -1, 1, -1, f"{180 - i * 10}")

    rt.gStyle.SetOptStat(0)
    # rt.gStyle.SetPalette(rt.kLightTemperature)

    canvas.SetLogx()
    hist.Draw("colz")

    input("Enter any symbol to quit: ")

    return


if __name__ == '__main__':
    plot_transmission_graph()
