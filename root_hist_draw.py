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


def create_hist_from_array(x: np.ndarray, y: np.ndarray, title: str) -> rt.TH1F:
    n = y.size
    hist = rt.TH1F(title, title, n - 1, x)

    for i in range(n):
        hist.Fill(x[i], y[i])

    hist.SetTitleOffset(0.2)

    size = .04

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetXaxis().SetTitleOffset(1.3)
    hist.GetXaxis().SetTitleSize(size)
    hist.GetXaxis().SetLabelSize(size)

    hist.GetYaxis().SetTitle("events per year per bin, 10^{-3}")
    hist.GetYaxis().SetTitleOffset(1.2)
    hist.GetYaxis().SetTitleSize(size)
    hist.GetYaxis().SetLabelSize(size)

    return hist


def draw_root_hist(sources: list[Source], source_numbers: list[int], energy: np.ndarray,
                   d_simple_reg: list[np.ndarray], simple_reg: list[np.ndarray], reg: list[np.ndarray]):
    colors = [602, 633, 419, 0]
    fill = [3654, 3645, 3095, 3001]

    canvas = rt.TCanvas("c", "c", 800, 600)
    canvas.SetLeftMargin(.13)
    canvas.SetBottomMargin(.13)
    canvas.SetRightMargin(.05)
    # canvas.SetTopMargin(.13)

    hist_simple, hist_complex, hist_d_simple = [], [], []
    title = ""

    for i, n in enumerate(source_numbers):
        source = sources[n]

        hist_simple.append(create_hist_from_array(energy, simple_reg[i] * 1e4, sources[n].name))
        hist_complex.append(create_hist_from_array(energy, reg[i] * 1e4, "1"))
        # hist_d_simple = create_hist_from_array(energy, d_simple_reg[i] * 1e4, sources[n].name)

        simple_integral = np.round(np.sum(hist_simple[i]) * 1e-3, 2)
        complex_integral = np.round(np.sum(hist_complex[i]) * 1e-3, 3)
        d_simple_integral = np.round(np.sum(hist_d_simple) * 1e-3, 2)

        hist_complex[i].SetLineColor(colors[i])
        hist_complex[i].SetFillColor(colors[i])
        hist_complex[i].SetLineWidth(2)
        hist_simple[i].SetLineColor(colors[i+1])
        hist_simple[i].SetLineWidth(2)
        # hist_d_simple.SetLineColor(colors[i+2])
        # hist_d_simple.SetFillColor(colors[i+2])

        if len(source_numbers) > 1:
            title = sources[n].name[:-3]
        else:
            title = sources[n].name

        hist_complex[i].SetFillStyle(fill[i])
        hist_simple[i].SetFillStyle(fill[i+1])
        hist_simple[i].SetFillStyle(fill[i+1])
        # hist_d_simple.SetFillStyle(fill[i+2])

        text = add_text(align_left=True)
        x_0, y_0, dy = .2, .75, .08

        if i == 0:
            # hist_simple[i].Draw("hist")
            # hist_d_simple.Draw("same hist")
            # hist_simple[i].Draw("same hist")
            hist_complex[i].Draw("same hist")
            text.DrawLatexNDC(x_0, 0.82, '#delta = ' + str(np.round(source.declination, 2)) + "#circ")
            # + ', ' + ' #Gamma = ' + str(source.gamma))
        else:
            # hist_simple[i].Draw("same hist")
            hist_complex[i].Draw("same hist")

        rt.gStyle.SetTitleAlign(11)
        rt.gStyle.SetTitleX(.99)

        # text.SetTextColor(colors[i + 2])
        # text.DrawLatexNDC(x_0, y_0 - (i + 2) * dy,
        #                   'R_{' + str(i + 3) + '} = ' + str(d_simple_integral) + ' #frac{counts}{year}')

        # text.SetTextColor(colors[i + 1])
        # text.DrawLatexNDC(x_0, y_0 - (i+1) * dy,
        #                   'R_{' + str(i + 2) + '} = ' + str(simple_integral) + ' #frac{counts}{year}')
        text.SetTextColor(colors[i])
        text.DrawLatexNDC(x_0, y_0 - i * dy,
                          'R_{' + str(i + 1) + '} = ' + str(complex_integral) + ' #frac{counts}{year}')

    text = add_text(24)
    text.DrawLatexNDC(.78, .945, 'Baikal-GVD MC, 10 clusters')
    text = add_text(26)
    text.DrawLatexNDC(.2, .95, title)

    rt.gStyle.SetOptStat(0)
    canvas.SetLogx()
    canvas.Update()

    input("Insert any symbol to quit: ")

    return


if __name__ == '__main__':
    print("Not for direct use!")
