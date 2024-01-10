import ROOT as rt
import matplotlib.pyplot as plt
import numpy as np


def check_root_transmission_table():

    file = rt.TFile.Open("data/data_mod.root", "READ")

    hist = file.Get("With_regeneration")

    xN, yN, zN = hist.GetNbinsX(), hist.GetNbinsY(), hist.GetNbinsZ()

    x1 = 100

    result = np.zeros([yN, zN])

    for i in range(yN):
        for j in range(zN):
            result[i, j] = hist.GetBinContent(x1 + 1, i + 1, j + 1)

    plt.imshow(result)
    plt.show()

    return


check_root_transmission_table()
