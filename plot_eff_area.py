import numpy as np
import ROOT as rt

from telescope import get_simple_telescope, get_Baikal


def draw_eff_area():

    t = get_Baikal("data/eff_area_trigger")

    m, n = t.angles.size, t.lg_energy.size

    a = np.linspace(-np.pi/2, 0, m)
    lg_e, e = t.lg_energy, t.energy

    xv = np.array(np.meshgrid(a, lg_e, indexing='ij')).T
    f_xv = t.ef_area(xv).T

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.1)
    canvas.SetBottomMargin(.1)
    canvas.SetRightMargin(.18)

    hist = rt.TH2F("Title", "Effective area on angle and energy", n-1, e, m-1, np.rad2deg(a))

    for i, a_i in enumerate(a):
        for j, e_j in enumerate(e):
            hist.Fill(e_j, np.rad2deg(a_i), f_xv[i, j] + 3e-6)

    size = .04

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetXaxis().SetTitleSize(size)
    hist.GetXaxis().SetLabelSize(size)
    hist.GetXaxis().SetTitleOffset(1.2)

    hist.GetYaxis().SetTitle("#theta, deg")
    hist.GetYaxis().SetTitleSize(size)
    hist.GetYaxis().SetLabelSize(size)
    hist.GetYaxis().SetTitleOffset(1.2)

    hist.GetZaxis().SetTitle("A_{ef}, m^{2}")
    hist.GetZaxis().SetTitleSize(size)
    hist.GetZaxis().SetLabelSize(size)
    hist.GetZaxis().SetTitleOffset(1.2)

    rt.gStyle.SetOptStat(0)
    hist.Draw("colz")
    canvas.SetLogx()
    canvas.SetLogz()

    input("Enter any symbol to quit: ")

    # plt.matshow(np.log10(f_xv.T + 1e-9), cmap='inferno')
    # q, p = 6, 3
    # plt.xticks(ticks=np.arange(0, lg_e.size, q), labels=np.round(lg_e[::q], 1))
    # plt.yticks(ticks=np.arange(0, a.size, p), labels=np.round(a[::p] / np.pi, 1))
    # plt.show()

    return


if __name__ == '__main__':
    draw_eff_area()
