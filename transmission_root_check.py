import ROOT as rt
import matplotlib.pyplot as plt
import numpy as np


def check_root_transmission_table():

    print("0 -> no regeneration, 1 -> with regeneration, 2 -> tau with regeneration")
    s = int(input("Enter a number 0, 1 or 2: "))
    
    file = rt.TFile.Open("data/data_mod.root", "READ")

    names = ["No_regeneration", "With_regeneration", "tau_with_regeneration"]

    if s != 0 and s != 1 and s != 2:
        print("Wrong value!")
        return

    hist = file.Get(names[s])

    xN, yN, zN = hist.GetNbinsX(), hist.GetNbinsY(), hist.GetNbinsZ()

    print("choose a zenith angle: 0 -> π/2, 179 -> π")
    x1 = int(input("Enter a number from 0 to 180: "))

    if x1 < 0 or x1 > 180:
        print("Wrong value!")
        return

    result = np.zeros([yN, zN])

    for i in range(yN):
        for j in range(zN):
            result[i, j] = hist.GetBinContent(x1 + 1, i + 1, j + 1)

    plt.title(f"transmission for z = {np.round(np.pi/2 *(1 + (x1+1)/180), 1)} rad")

    plt.xlabel(r"$E_{in},~GeV$")
    plt.ylabel(r"$E_{out},~GeV$")

    plt.xticks(np.linspace(0, 198, 8), np.round(np.linspace(3, 10, 8), 0))
    plt.yticks(np.linspace(0, 198, 8), np.round(np.linspace(3, 10, 8), 0))

    plt.tight_layout()

    plt.imshow(result, origin='lower')
    plt.show()

    return


check_root_transmission_table()
