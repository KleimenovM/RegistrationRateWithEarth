import matplotlib.pyplot as plt
import numpy as np
import ROOT as rt

from tools import smart_division
from transmission_function import TransmissionFunction


def plot_transmission_graph(gamma=2):
    tf = TransmissionFunction()

    angles = np.linspace(-np.pi / 2, 0, 180)
    e = tf.energy

    flux = e**(-gamma)

    m, n = angles.size, e.size

    hist = rt.TH2F("Title", "Transmission matrix for #Gamma = 2", n-1, e, m-1, angles)

    for i, a in enumerate(angles):
        f_i = smart_division(tf.convolution(a, flux, 2), flux)
        for j, e_j in enumerate(e):
            hist.Fill(e_j, a, f_i[j])

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.1)
    canvas.SetBottomMargin(.1)
    canvas.SetRightMargin(.18)

    size = 0.04

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetXaxis().SetTitleSize(size)
    hist.GetXaxis().SetLabelSize(size)
    hist.GetXaxis().SetTitleOffset(1.2)

    hist.GetYaxis().SetTitle("#theta, rad")
    hist.GetYaxis().SetTitleSize(size)
    hist.GetYaxis().SetLabelSize(size)
    hist.GetYaxis().SetTitleOffset(1.2)

    hist.GetZaxis().SetTitle("#phi/#phi_{0}")
    hist.GetZaxis().SetTitleSize(size)
    hist.GetZaxis().SetLabelSize(size)
    hist.GetZaxis().SetTitleOffset(1.2)

    rt.gStyle.SetOptStat(0)
    # rt.gStyle.SetPalette(rt.kLightTemperature)

    canvas.SetLogx()
    hist.Draw("colz")

    input("Enter any symbol to quit: ")

    return


if __name__ == '__main__':
    plot_transmission_graph()
