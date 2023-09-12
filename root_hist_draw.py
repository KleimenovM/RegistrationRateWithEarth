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
    print(f"n={y.size}, n+1={x.size}")
    print(x)
    hist = rt.TH1F(title, title, n, x)

    for i in range(n):
        hist.Fill(x[i], y[i])

    hist.GetXaxis().SetTitle("E, GeV")
    hist.GetYaxis().SetTitle("counts per year")

    return hist


def draw_root_hist(sources: list[Source], source_numbers: list[int],
                   energy: np.ndarray, simple_reg: list[np.ndarray], reg: list[np.ndarray]):
    colors = [602, 633, 15]
    fill = [3354, 3345, 3001]

    canvas = rt.TCanvas()

    for i, n in enumerate(source_numbers):
        source = sources[n]
        hist_simple = create_hist_from_array(energy, simple_reg[i], "1")
        hist_complex = create_hist_from_array(energy, reg[i], "2")

        simple_integral = np.round(np.sum(hist_simple) * 14, 3)
        complex_integral = np.round(np.sum(hist_complex) * 14, 3)

        hist_complex.SetFillColor(colors[i])
        hist_simple.SetLineColor(colors[i+1])

        if len(source_numbers) > 1:
            hist_complex.SetTitle(sources[n].name[:-3])
        else:
            hist_complex.SetTitle(sources[n].name)

        if i < 2:
            hist_complex.SetFillStyle(fill[i])
            hist_simple.SetFillStyle(fill[i+1])

        text = add_text()

        if i == 0:
            hist_complex.Draw("hist")
            hist_simple.Draw("same hist")
            text.DrawLatexNDC(0.25, 0.85, '#delta=' + str(np.round(source.declination, 2)))
        else:
            hist_complex.Draw("same hist")

        text.SetTextColor(colors[i])
        text.DrawLatexNDC(0.25, 0.75 - i * .1, 'R_{' + str(i + 1) + '} = ' + str(complex_integral) + ' #frac{counts}{year}')
        text.SetTextColor(colors[i+1])
        text.DrawLatexNDC(0.25, 0.75 - (i + 1) * .1, 'R_{' + str(i + 2) + '} = ' + str(simple_integral) + ' #frac{counts}{year}')
    rt.gStyle.SetOptStat(0)
    canvas.SetLogx()
    canvas.Update()

    h = input("Insert any symbol to quit: ")

    return


if __name__ == '__main__':
    print("Not for direct use!")
