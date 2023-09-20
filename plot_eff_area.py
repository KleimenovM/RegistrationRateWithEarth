import numpy as np
import matplotlib.pyplot as plt
import ROOT as rt

from main import get_Baikal
from telescope import get_simple_telescope


def draw_eff_area():

    t = get_Baikal()
    t2 = get_simple_telescope("BaikalGVD", [51, 46], "data/eff_area.root")

    m, n = t.angles.size, t.lg_energy.size

    a = np.linspace(-np.pi/2, 0, m)
    print(a)
    lg_e, e = t.lg_energy, t.energy

    xv = np.array(np.meshgrid(a, lg_e, indexing='ij')).T
    f_xv = t.ef_area(xv).T

    canvas = rt.TCanvas()
    hist = rt.TH2F("Title", "Effective area on angle and energy", n-1, e, m, np.linspace(-11*np.pi/20, 0, m+1))

    for i, a_i in enumerate(a):
        for j, e_j in enumerate(e):
            hist.Fill(e_j, a_i - np.pi/20, f_xv[i, j])

    yv = np.array(np.meshgrid(-np.pi/2, lg_e, indexing='ij')).T
    f_yv = t2.ef_area(yv).T

    hist2 = rt.TH2F("Title2", "Title2", n - 1, e, m, np.linspace(-11 * np.pi / 20, 0, m + 1))

    for j, e_j in enumerate(e):
        hist2.Fill(e_j, -np.pi/40, f_yv[0, j])

    hist2.SetFillColor(rt.kRed)

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetYaxis().SetTitle("#theta, rad")

    rt.gStyle.SetOptStat(0)
    hist.Draw("LEGO2")
    hist2.Draw("same LEGO")
    canvas.SetLogx()
    # canvas.SetLogz()

    input("Enter any symbol to quit: ")

    # plt.matshow(np.log10(f_xv.T + 1e-9), cmap='inferno')
    # q, p = 6, 3
    # plt.xticks(ticks=np.arange(0, lg_e.size, q), labels=np.round(lg_e[::q], 1))
    # plt.yticks(ticks=np.arange(0, a.size, p), labels=np.round(a[::p] / np.pi, 1))
    # plt.show()

    return


if __name__ == '__main__':
    draw_eff_area()
