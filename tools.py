# Neutrino through Earth propagation
# Tools description

import numpy as np
from scipy.interpolate import interp1d


def deg_to_rad(deg: list):
    result = .0
    for i, a in enumerate(deg):
        result += deg[i] / (60 ** i) / 180 * np.pi
    return result


def rot_matrix(rotation_angle: float):
    cp, sp = np.cos(rotation_angle), np.sin(rotation_angle)
    return np.array([[1, 0, 0], [0, sp, -cp], [0, cp, sp]])


def sph_coord(r, theta, phi):
    x = r * np.cos(theta) * np.cos(phi)
    y = r * np.cos(theta) * np.sin(phi)
    z = r * np.sin(theta)
    return x, y, z


def extrapolating_spline(x, x0, y0, if_delta=False):
    if if_delta:
        spline = interp1d(np.log10(x0), y0, kind='linear', fill_value='extrapolate')
        return spline(np.log10(x))
    spline = interp1d(np.log10(x0), np.log(y0), kind='linear', fill_value='extrapolate')
    y_p = spline(np.log10(x))
    # y_p[y_p > 0] = 0
    return np.exp(y_p)


if __name__ == '__main__':
    "Not for direct use"
