# Neutrino through Earth propagation
# Tools description

import numpy as np
from scipy.interpolate import interp1d

Ag = np.array([[-0.0548755601367195, -0.8734370902532698, -0.4838350155472244],
               [+0.4941094280132430, -0.4448296298016944, +0.7469822445004389],
               [-0.8676661489582886, -0.1980763737056720, +0.4559837761713720]])


def deg_to_rad(deg: list):
    result = .0
    for i in range(len(deg)):
        result += deg[i] / (60 ** i) / 180 * np.pi
    return result


def hours_to_rad(h):
    return h / 12 * np.pi


def rad_to_hours(r):
    return r / np.pi * 12


def get_a_const_func(value):
    def c_func(*args):
        return value

    return c_func()


def rot_matrix(rotation_angle: float):
    cp, sp = np.cos(rotation_angle), np.sin(rotation_angle)
    return np.array([[1, 0, 0],
                     [0, sp, -cp],
                     [0, cp, sp]])


def sph_coord(r, theta, phi):
    x = r * np.cos(theta) * np.cos(phi)
    y = r * np.cos(theta) * np.sin(phi)
    z = r * np.sin(theta)
    return x, y, z


def galactic2equatorial(ll, b, get_delta=False, get_alpha=False):
    r_gal = sph_coord(1, theta=b, phi=ll)
    r_eq = Ag.T.dot(r_gal)

    # get all angles
    if get_delta and get_alpha:
        delta = np.arcsin(r_eq[2])
        c_delta = np.cos(delta)
        c_alpha = r_eq[0] / c_delta
        s_alpha = r_eq[1] / c_delta
        if s_alpha >= 0:
            alpha = np.arccos(c_alpha)
        else:
            alpha = 2 * np.pi - np.arccos(c_alpha)
        return delta, alpha  # delta (rad), ra (rad)

    # get zenith angle only
    if get_delta:
        return np.arcsin(r_eq[2])  # delta (rad)

    # 3D-coordinates
    return r_eq


def equatorial2galactic(delta, alpha):
    r_eq = sph_coord(1, theta=delta, phi=alpha)
    r_gal = Ag.dot(r_eq)
    return r_gal


def smart_division(a, b):
    good_indices = (b != 0.)
    result = np.zeros(b.size)
    result[good_indices] = a[good_indices] / b[good_indices]
    return result


def extrapolating_spline(x, x0, y0, if_delta=False):
    if if_delta:
        spline = interp1d(np.log10(x0), y0, kind='linear', fill_value='extrapolate')
        return spline(np.log10(x))
    spline = interp1d(np.log10(x0), np.log(y0), kind='linear', fill_value='extrapolate')
    y_p = spline(np.log10(x))
    # y_p[y_p > 0] = 0
    return np.exp(y_p)


if __name__ == '__main__':
    print("Not for direct use")
