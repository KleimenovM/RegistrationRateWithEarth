import numpy as np
import matplotlib.pyplot as plt


from source import get_sources


def coordinates_conversion(delta, phi):
    # See https://en.wikipedia.org/wiki/Hammer_projection
    delta = delta / 24 * 2 * np.pi
    phi = np.deg2rad(phi)
    denominator = np.sqrt(1 + np.cos(phi) * np.cos(delta/2))
    x = 2**3/2 * np.cos(phi) * np.sin(delta/2) / denominator
    y = 2**0.5 * np.sin(phi) / denominator
    return x, y


def plot_grid(delta_step=30, phi_step=4):
    phi_brd = [-90, 90]
    lambda_brd = [-12, 12]

    # draw parallels
    d = phi_brd[0]
    ph = np.linspace(lambda_brd[0], lambda_brd[1], 100)
    while d <= phi_brd[1]:
        x, y = coordinates_conversion(ph, d)
        plt.plot(x, y, color='gray', alpha=.5)

        # text
        x_t, y_t = coordinates_conversion(-12, d)
        if d < 0:
            y_t -= .14
            x_t -= .2
        if d == 0:
            y_t -= .07
            x_t -= .1
        plt.text(x_t - 0.2, y_t + .02, str(d) + r'$^\circ$', fontsize=14)

        d += delta_step

    # draw meridians
    d = lambda_brd[0]
    phi = np.linspace(phi_brd[0], phi_brd[1], 100)
    while d <= lambda_brd[1]:
        x, y = coordinates_conversion(d, phi)
        plt.plot(x, y, color='gray', alpha=.5)

        # text
        x_t, y_t = coordinates_conversion(d, 0)
        plt.text(x_t + .05, y_t - .14, str(12 - d) + r'$^h$', fontsize=14)

        d += phi_step

    x1, y1 = coordinates_conversion(-12, phi)
    x2, y2 = coordinates_conversion(12, phi)
    plt.fill_betweenx(y2, x1, x2, color='#009', alpha=.1)

    print(x1, x2)

    return


def draw_hammer_aitoff():

    plt.figure(figsize=(10, 6))

    plot_grid()
    sources = get_sources("data/source_table.csv")

    for i, s in enumerate(sources):
        if_change_name = 0
        if i > 0 and s.name[:-3] == sources[i-1].name[:-3]:
            continue
        if i < len(sources) - 1 and s.name[:-3] == sources[i+1].name[:-3]:
            name = s.name[:-3]
        else:
            name = s.name
        x, y = coordinates_conversion(12 - s.right_ascension, s.declination)
        print(x, y)
        plt.scatter(x, y, label=name)
        if s.name == "Vela Jr":
            plt.text(x - 0.2, y - 0.15, name, fontsize=12)
        else:
            plt.text(x - 0.2, y + 0.05, name, fontsize=12)

    plt.axis('off')
    plt.tight_layout()
    plt.show()

    return


if __name__ == '__main__':
    draw_hammer_aitoff()
