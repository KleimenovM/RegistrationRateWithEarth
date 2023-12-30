import ROOT as rt
import numpy as np

from source import Source


def add_text(fs: int = 22, align_left=False):
    text = rt.TLatex(.5, .5, "")
    text.SetTextFont(43)
    text.SetTextSize(fs)
    if align_left:
        text.SetTextAlign(11)
    else:
        text.SetTextAlign(22)
    return text


def create_hist_from_array(x: np.ndarray, y: np.ndarray, title: str, x0=None) -> rt.TH1F:
    if x0 is None:
        x0 = x

    m = x0.size

    n = y.size
    hist = rt.TH1F(title, title, m - 1, x0)

    for i in range(n):
        hist.Fill(x[i], y[i])

    hist.SetTitleOffset(0.2)

    size = .04

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetXaxis().SetTitleOffset(1.5)
    hist.GetXaxis().SetLabelOffset(0.015)
    hist.GetXaxis().SetTitleSize(size)
    hist.GetXaxis().SetLabelSize(size)

    hist.GetYaxis().SetTitle("events per 5 years per bin")
    hist.GetYaxis().SetTitleOffset(1.5)
    hist.GetYaxis().SetTitleSize(size)
    hist.GetYaxis().SetLabelSize(size)

    return hist


def draw_root_hist(sources: list[Source], source_numbers: list[int],
                   energy_s: np.ndarray, energy_c: np.ndarray = None,
                   simple_reg: list[np.ndarray] = None, reg: list[np.ndarray] = None,
                   value_c: float = 1, value_s: float = 1, caption_pos='left'):
    energy_s = energy_s[1:]
    if energy_c is None:
        energy_c = energy_s
    else:
        energy_c = energy_c[1:]

    colors = [602, 633, 419, 0]
    fill = [3654, 3645, 3095, 3001]

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.13)
    canvas.SetBottomMargin(.13)
    canvas.SetRightMargin(.13)
    canvas.SetTopMargin(.13)

    hist_simple, hist_complex, hist_d_simple = [], [], []
    title = ""

    for i, n in enumerate(source_numbers):
        source = sources[n]

        hist_simple.append(create_hist_from_array(energy_s, simple_reg[i][1:] * value_s, sources[n].name, x0=energy_c))
        hist_complex.append(create_hist_from_array(energy_c, reg[i][1:] * value_c, "1"))

        simple_integral = np.round(np.sum(hist_simple[i]), 2)
        complex_integral = np.round(np.sum(hist_complex[i]), 2)

        hist_complex[i].SetLineColor(colors[i])
        hist_complex[i].SetFillColor(colors[i])
        hist_complex[i].SetLineWidth(2)
        hist_simple[i].SetLineColor(colors[i + 1])
        hist_simple[i].SetLineWidth(2)

        if len(source_numbers) > 1:
            title = sources[n].name[:-3]
        else:
            title = sources[n].name

        hist_complex[i].SetFillStyle(fill[i])
        hist_simple[i].SetFillStyle(fill[i + 1])
        hist_simple[i].SetFillStyle(fill[i + 1])

        text = add_text(align_left=True)

        if caption_pos == 'left':
            x_0 = .17
        else:
            x_0 = .7

        y_0, dy = .75, .08

        if i == 0:
            if hist_simple[i].GetMaximum() > hist_complex[i].GetMaximum():
                hist_simple[i].Draw("hist")
                hist_complex[i].Draw("same hist")
            else:
                hist_complex[i].Draw("hist")
                hist_simple[i].Draw("same hist")
            text.DrawLatexNDC(x_0, y_0 + 0.07, '#delta = ' + str(np.round(source.declination, 2)) + "#circ")
        else:
            hist_simple[i].Draw("same hist")
            hist_complex[i].Draw("same hist")

        ratio = hist_complex[i] / hist_simple[i]
        rightmax = ratio.GetMaximum() * 1.05

        axis = rt.TGaxis(energy_c.max(), 0.0,
                         energy_c.max(), hist_simple[i].GetMaximum() * 1.05,
                         0.0, rightmax, 12, "+L")
        axis.SetTitle("reconstruction / trigger")
        # axis.SetLabelSize(0.03)
        # axis.SetTitleSize(0.03)
        axis.SetLabelFont(42)
        axis.SetTitleFont(42)
        axis.SetTitleOffset(1.3)

        axis.Draw()

        scale = hist_simple[i].GetMaximum() / rightmax
        ratio.SetLineColor(rt.kBlack)
        ratio.SetFillColor(rt.kWhite)
        ratio.SetFillStyle(3001)
        ratio.SetLineWidth(3)
        ratio.Scale(scale)
        ratio.Draw("same hist")
        axis.Draw()

        # hist_simple[i].GetXaxis().SetRangeUser(1e2, 1e7)
        # hist_complex[i].GetXaxis().SetRangeUser(1e2, 1e7)

        rt.gStyle.SetTitleAlign(11)
        rt.gStyle.SetTitleX(.99)

        text.SetTextColor(colors[i + 1])
        text.DrawLatexNDC(x_0, y_0 - (i + 1) * dy,
                          'R_{trig} = ' + str(simple_integral) + ' #frac{counts}{5 years}')
        text.SetTextColor(colors[i])
        text.DrawLatexNDC(x_0, y_0 - i * dy,
                          'R_{reco} = ' + str(complex_integral) + ' #frac{counts}{5 years}')

    text = add_text(24, align_left=True)
    text.SetTextColor(colors[0])
    # text.DrawLatexNDC(.5, .955, 'Baikal-GVD MC, 20 clusters, 5 yr')
    text.DrawLatexNDC(.3, .955, 'Baikal-GVD MC, 20 clusters, 5 yr, reconstruction')
    text.SetTextColor(colors[1])
    text.DrawLatexNDC(.3, .905, 'Baikal-GVD MC, 20 clusters, 5 yr, trigger')
    # text.DrawLatexNDC(.5, .905, 'KM3Net, Full Telescope, 2 blocks, 5 yr')
    text = add_text(26)
    text.DrawLatexNDC(.15, .95, title)

    rt.gStyle.SetOptStat(0)
    canvas.SetLogx()
    canvas.SetGrayscale()
    canvas.Update()

    input("Insert any symbol to quit: ")

    return


def draw_root_ext(energy_s: np.ndarray, energy_c: np.ndarray = None,
                  simple_reg: list[np.ndarray] = None, reg: list[np.ndarray] = None,
                  value_c: float = 1, value_s: float = 1, caption_pos='left'):
    energy_s = energy_s[1:]
    if energy_c is None:
        energy_c = energy_s
    else:
        energy_c = energy_c[1:]

    colors = [602, 633, 419, 0]
    fill = [3654, 3645, 3095, 3001]

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.13)
    canvas.SetBottomMargin(.13)
    canvas.SetRightMargin(.05)
    canvas.SetTopMargin(.13)

    title = "Galactic ridge"

    hist_simple = create_hist_from_array(energy_s, simple_reg[1:] * value_s, "")  #, x0=energy_c)
    hist_complex = create_hist_from_array(energy_c, reg[1:] * value_c, "1")

    simple_integral = np.round(np.sum(hist_simple), 2)
    complex_integral = np.round(np.sum(hist_complex), 2)

    hist_complex.SetLineColor(colors[0])
    hist_complex.SetFillColor(colors[0])
    hist_complex.SetLineWidth(2)
    hist_simple.SetLineColor(colors[1])
    hist_simple.SetLineWidth(2)

    hist_complex.SetFillStyle(fill[0])
    hist_simple.SetFillStyle(fill[1])

    text = add_text(align_left=True)

    if caption_pos == 'left':
        x_0 = .2
    else:
        x_0 = .67

    y_0, dy = .77, .08

    if hist_simple.GetMaximum() > hist_complex.GetMaximum():
        hist_simple.Draw("hist")
        hist_complex.Draw("same hist")
    else:
        hist_complex.Draw("hist")
        hist_simple.Draw("same hist")

    rt.gStyle.SetTitleAlign(11)
    rt.gStyle.SetTitleX(.99)

    text.SetTextColor(colors[1])
    text.DrawLatexNDC(x_0, y_0 - dy, 'R_{trig} = ' + str(simple_integral) + ' #frac{counts}{5 years}')
    text.SetTextColor(colors[0])
    text.DrawLatexNDC(x_0, y_0, 'R_{reco} = ' + str(complex_integral) + ' #frac{counts}{5 years}')

    text = add_text(24, align_left=True)
    text.SetTextColor(colors[0])
    text.DrawLatexNDC(.5, .955, 'Baikal-GVD MC, 20 clusters, 5 yr')
    # text.DrawLatexNDC(.3, .955, 'Baikal-GVD MC, 20 x 1 clusters, 5 yr, reconstruction')
    text.SetTextColor(colors[1])
    # text.DrawLatexNDC(.3, .905, 'Baikal-GVD MC, 4 x 5 clusters, 5 yr, trigger')
    text.DrawLatexNDC(.5, .905, 'KM3Net, Full Telescope, 2 blocks, 5 yr')
    text = add_text(26)
    text.DrawLatexNDC(.15, .95, title)

    rt.gStyle.SetOptStat(0)
    canvas.SetLogx()
    # canvas.SetGrayscale()
    canvas.Update()

    input("Insert any symbol to quit: ")

    return


if __name__ == '__main__':
    print("Not for direct use!")
