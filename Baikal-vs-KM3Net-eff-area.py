import numpy as np
import matplotlib.pyplot as plt
import ROOT as rt

from telescope import Telescope, get_simple_telescope_from_txt, get_simple_telescope
from tools import deg_to_rad


def pyplot_figure(energy: np.ndarray,
                  area1: np.ndarray, area2: np.ndarray,
                  name1: str, name2: str):
    fig, ax1 = plt.subplots(figsize=(8, 6))

    ax1.plot(energy, area1, color='blue', label=name1)
    ax1.plot(energy, area2, color='red', label=name2)

    plt.xlim(1e2, 1e6)
    plt.ylim(1e-3, 1e3)
    plt.xscale('log')
    plt.yscale('log')

    ax2 = ax1.twinx()
    ax2.plot(energy, area1 / area2, color='black')
    ax2.plot(energy, np.ones(energy.size), color='gray', linestyle='dashed')

    plt.xlim(1e2, 1e6)
    plt.ylim(1e-1, 1e1)
    plt.xscale('log')
    plt.yscale('log')

    plt.show()
    return


def root_figure(energy: np.ndarray,
                area1: np.ndarray, area2: np.ndarray,
                name1: str, name2: str):

    n = energy.size
    title = "Baikal-GVD & KM3Net trigger effective area"
    hist1 = rt.TH1F(name1, title, n-1, energy)
    hist2 = rt.TH1F(name2, title, n-1, energy)
    hist3 = rt.TH1F("ratio", title, n-1, energy)
    hist4 = rt.TH1F("unit", title, n-1, energy)

    for i, e_i in enumerate(energy):
        hist1.Fill(e_i, area1[i])
        hist2.Fill(e_i, area2[i])
        hist3.Fill(e_i, area1[i] / area2[i])
        hist4.Fill(e_i, 1.0)

    colors = [rt.kBlue, rt.kRed, 419, 0]
    fill = [3654, 3945, 3095, 3001]

    canvas = rt.TCanvas("c", "c", 800, 800)
    canvas.Draw()
    # canvas.SetLeftMargin(.12)
    # canvas.SetBottomMargin(.1)
    # canvas.SetRightMargin(.05)
    canvas.SetTopMargin(.05)

    pad1 = rt.TPad("p1", "p1", 0, 0.4, 1, 1)
    pad1.Draw()
    legend = rt.TLegend(0.15, 0.75, 0.6, 0.88)

    for i, hist in enumerate([hist1, hist2]):
        size = .045

        hist.GetXaxis().SetTitle("E, GeV")
        # hist.GetXaxis().SetTitleOffset(1.0)
        hist.GetXaxis().SetTitleSize(size)
        hist.GetXaxis().SetLabelSize(size)
        hist.GetXaxis().SetLabelOffset(0.03)

        hist.GetYaxis().SetTitle("effective area, m^{2}")
        hist.GetYaxis().SetTitleOffset(1)
        hist.GetYaxis().SetTitleSize(size)
        hist.GetYaxis().SetLabelSize(size)

        legend.AddEntry(hist, hist.GetName())

        hist.SetLineWidth(2)
        hist.SetLineColor(colors[i])
        hist.SetFillColor(colors[i])
        hist.SetFillStyle(fill[i])

    rt.gStyle.SetOptStat(0)
    pad1.SetLogx()
    pad1.SetLogy()
    pad1.SetBottomMargin(0)

    pad1.cd()
    hist1.Draw("hist")
    hist2.Draw("same hist")
    legend.Draw()

    rt.gPad.SetTicky(2)
    rt.gPad.SetTickx(2)

    # PAD 2
    canvas.cd()
    pad2 = rt.TPad("p2", "p2", 0, 0, 1, 0.4)
    pad2.Draw()
    rt.gStyle.SetOptTitle(0)

    pad2.SetLogx()
    hist3.SetLineWidth(2)
    hist3.SetLineColor(colors[2])

    hist4.SetLineWidth(2)
    hist4.SetLineStyle(9)
    hist4.SetLineColor(rt.kBlack)

    size = 0.07
    hist3.GetXaxis().SetTitle("E, GeV")
    hist3.GetXaxis().SetTitleOffset(0.5)
    hist3.GetXaxis().SetTitleSize(size)
    hist3.GetXaxis().SetLabelSize(size)

    hist3.GetYaxis().SetTitle("A_{Baikal} / A_{KM3NeT}")
    hist3.GetYaxis().SetTitleOffset(0.6)
    hist3.GetYaxis().SetTitleSize(size)
    hist3.GetYaxis().SetLabelSize(size)

    pad2.cd()

    line = rt.TLine(100, 1, 800, 1)
    line.Draw()
    rt.gPad.SetTicky(2)
    # rt.gPad.SetTickx(2)

    pad2.SetTopMargin(.0)
    pad2.SetBottomMargin(1.2)
    hist3.Draw("hist")
    hist4.Draw("same hist")

    # canvas.SetGrayscale()
    canvas.Update()

    input("Type anything to exit: ")

    return


def compare_eff_areas():
    """
    Trigger level effective area comparison
    @return: 0
    """

    # Baikal
    baikal: Telescope = get_simple_telescope(name="Baikal-GVD-trigger", latitude=[51, 46],
                                             filename="data/eff_area_single_cluster.root",
                                             histname="hnu_trigger",
                                             brd_angle=[30])

    baikal_clusters = 20  # 20 Baikal clusters

    # KM3Net
    km3net: Telescope = get_simple_telescope_from_txt(name="KM3Net-trigger", latitude=[36, 16],
                                                      filename="data/KM3Net-total.txt",
                                                      brd_angle=[-12])

    b_lg_e, k_lg_e = baikal.lg_energy, km3net.lg_energy
    lg_e_min, lg_e_max = max(min(b_lg_e), min(k_lg_e)), min(max(b_lg_e), max(k_lg_e), 5.9)
    lg_energy = np.arange(lg_e_min, lg_e_max, b_lg_e[1] - b_lg_e[0])
    energy = 10**lg_energy
    theta = deg_to_rad([-90])

    grid_x, grid_y = np.meshgrid(theta, lg_energy, indexing='ij')
    grid = np.array([grid_x, grid_y]).T

    baikal_ef_area = baikal.ef_area(grid).T[0] * baikal_clusters
    km3net_ef_area = km3net.ef_area(grid).T[0]

    name1, name2 = "Baikal-GVD (20 clusters, trigger)", "KM3Net (2 blocks, trigger)"

    # Matplotlib realization
    # pyplot_figure(energy=energy,
    #               area1=baikal_ef_area, name1=name1,
    #               area2=km3net_ef_area, name2=name2)

    # ROOT realization
    root_figure(energy=energy,
                area1=baikal_ef_area, name1=name1,
                area2=km3net_ef_area, name2=name2)

    return 0


if __name__ == '__main__':
    compare_eff_areas()
