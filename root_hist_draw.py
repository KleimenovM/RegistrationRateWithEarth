import ROOT as rt
import numpy as np

from source import Source


def add_text():
    text = rt.TLatex(.5, .5, "")
    text.SetTextAlign(22)
    text.SetTextFont(43)
    text.SetTextSize(20)
    return text


def create_hist_from_array(x: np.ndarray, y: np.ndarray, title: str) -> rt.TH1F:
    n = y.size
    hist = rt.TH1F(title, title, n, x)

    for i in range(n):
        hist.Fill(x[i], y[i])

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetYaxis().SetTitle("counts per year")

    return hist


def draw_root_hist(sources: list[Source], source_numbers: list[int], energy: np.ndarray,
                   d_simple_reg: list[np.ndarray], simple_reg: list[np.ndarray], reg: list[np.ndarray]):
    colors = [602, 633, 419]
    fill = [3654, 3355, 3645]

    canvas = rt.TCanvas()

    for i, n in enumerate(source_numbers):
        source = sources[n]
        hist_simple = create_hist_from_array(energy, simple_reg[i], sources[n].name)
        hist_complex = create_hist_from_array(energy, reg[i], "2")
        hist_d_simple = create_hist_from_array(energy, d_simple_reg[i], "3")

        simple_integral = np.round(np.sum(hist_simple) * 14, 3)
        complex_integral = np.round(np.sum(hist_complex) * 14, 3)
        d_simple_integral = np.round(np.sum(hist_d_simple) * 14, 3)

        hist_complex.SetFillColor(colors[i])
        hist_simple.SetLineColor(colors[i+1])
        hist_d_simple.SetLineColor(colors[i + 2])
        hist_d_simple.SetFillColor(colors[i+2])

        if len(source_numbers) > 1:
            hist_complex.SetTitle(sources[n].name[:-3])
        else:
            hist_complex.SetTitle(sources[n].name)

        if i < 2:
            hist_complex.SetFillStyle(fill[i])
            hist_simple.SetFillStyle(fill[i+1])
            hist_d_simple.SetFillStyle(fill[i+2])

        text = add_text()

        if i == 0:
            hist_simple.Draw("hist")
            hist_complex.Draw("same hist")
            hist_d_simple.Draw("same_hist")
            text.DrawLatexNDC(0.25, 0.85, '#delta = ' + str(np.round(source.declination, 2))
                              + ', ' + ' #Gamma = ' + str(source.gamma))
        else:
            hist_complex.Draw("same hist")

        text.SetTextColor(colors[i])
        text.DrawLatexNDC(0.25, 0.75 - i * .1, 'R_{' + str(i + 1) + '} = ' + str(complex_integral) + ' #frac{counts}{year}')
        text.SetTextColor(colors[i + 1])
        text.DrawLatexNDC(0.25, 0.75 - (i + 1) * .1, 'R_{' + str(i + 2) + '} = ' + str(simple_integral) + ' #frac{counts}{year}')
        text.SetTextColor(colors[i + 2])
        text.DrawLatexNDC(0.25, 0.75 - (i + 2) * .1,
                          'R_{' + str(i + 3) + '} = ' + str(d_simple_integral) + ' #frac{counts}{year}')
    rt.gStyle.SetOptStat(0)
    canvas.SetLogx()
    canvas.Update()

    input("Insert any symbol to quit: ")

    return


if __name__ == '__main__':
    print("Not for direct use!")
