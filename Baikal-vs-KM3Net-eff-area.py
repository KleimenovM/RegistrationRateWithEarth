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

    for i, e_i in enumerate(energy):
        hist1.Fill(e_i, area1[i])
        hist2.Fill(e_i, area2[i])

    colors = [602, 633, 419, 0]
    fill = [3654, 3645, 3095, 3001]

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.12)
    canvas.SetBottomMargin(.1)
    canvas.SetRightMargin(.05)
    canvas.SetTopMargin(.1)

    legend = rt.TLegend(0.15, 0.75, 0.6, 0.88)

    for i, hist in enumerate([hist1, hist2]):
        size = .04

        hist.GetXaxis().SetTitle("E, GeV")
        hist.GetXaxis().SetTitleOffset(1.0)
        hist.GetXaxis().SetTitleSize(size)
        hist.GetXaxis().SetLabelSize(size)
        hist.GetXaxis().SetLimits(10**2, 10**6)

        hist.GetYaxis().SetTitle("effective area, m^{2}")
        hist.GetYaxis().SetTitleOffset(1.2)
        hist.GetYaxis().SetTitleSize(size)
        hist.GetYaxis().SetLabelSize(size)

        legend.AddEntry(hist, hist.GetName())

        hist.SetLineWidth(2)
        hist.SetLineColor(colors[i])
        hist.SetFillColor(colors[i])
        hist.SetFillStyle(fill[i])

    rt.gStyle.SetOptStat(0)
    canvas.SetLogx()
    canvas.SetLogy()

    hist1.Draw("hist")
    hist2.Draw("same hist")
    legend.Draw()
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
    pyplot_figure(energy=energy,
                  area1=baikal_ef_area, name1=name1,
                  area2=km3net_ef_area, name2=name2)

    # ROOT realization
    root_figure(energy=energy,
                area1=baikal_ef_area, name1=name1,
                area2=km3net_ef_area, name2=name2)

    return 0


if __name__ == '__main__':
    compare_eff_areas()
