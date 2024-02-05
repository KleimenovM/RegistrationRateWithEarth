import numpy as np
from matplotlib import pyplot as plt

from source import get_sources
from source_extended import ExtendedSource, galactic_center
from tools import galactic2equatorial, rad_to_hours, hours_to_rad


def draw_hammer_ext():
    plt.figure(figsize=(10, 6))

    plt.subplot(111, projection="hammer")

    # Galactic plane
    m = 1000
    points_ll, points_b = np.linspace(0, 2 * np.pi, m), np.zeros(m)

    x, y = np.zeros(m), np.zeros(m)
    for i in range(m):
        y[i], x[i] = galactic2equatorial(points_ll[i], points_b[i], get_delta=True, get_alpha=True)

    plt.scatter(np.pi - x, y, 0.3, color='black', alpha=.8, label='Galactic plane')

    # Galactic center
    source: ExtendedSource = galactic_center(num=1000)

    d_c, ra_c = source.equatorial_dec_grid, source.equatorial_ra_grid
    n = source.num

    x_s, y_s = np.zeros(n), np.zeros(n)
    for i in range(n):
        x_s[i] = np.pi - ra_c[i]
        y_s[i] = d_c[i]
    plt.scatter(x_s, y_s, 3, color='orange', alpha=0.8, label='Galactic center')

    # depict the sources
    sources = get_sources("data/source_table.csv")

    for i, s in enumerate(sources):
        if i > 0 and s.name[:-3] == sources[i - 1].name[:-3]:
            continue
        if i < len(sources) - 1 and s.name[:-3] == sources[i + 1].name[:-3]:
            name = s.name[:-3]
        else:
            name = s.name
        x, y = np.pi - hours_to_rad(s.right_ascension), s.delta

        plt.scatter(x, y, label=name)
        print(s.name)
        if s.name == "Vela Jr" or s.name == "NGC 1068":
            plt.text(x - 0.05, y - 0.15, name, fontsize=12)
        elif s.name == "Vela X":
            plt.text(x + 0.1, y - 0.03, name, fontsize=12)
        else:
            plt.text(x + 0.02, y + 0.08, name, fontsize=12)

    plt.grid(True)
    plt.tight_layout()
    m = 6
    ticks = [str(int(t))+r"${}^h$" for t in np.linspace(24 * (m-1)/m, 24*1/m, m-1)]
    plt.xticks(np.linspace(-np.pi * (m-2)/m, np.pi * (m-2)/m, m-1), ticks,
               fontsize=12, verticalalignment='bottom')
    plt.yticks(fontsize=14)
    plt.show()

    return


if __name__ == '__main__':
    # draw_hammer_aitoff_sources()
    draw_hammer_ext()
