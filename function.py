import math
import numpy as np
from scipy.optimize import curve_fit
import scipy.interpolate as interpolate


def load_dipole_function(file="Dipgrid.txt"):
    list_pos_xy = []
    list_values1 = []
    list_values2 = []
    list_values3 = []
    with open(file, "r") as file:
        cur_deg = 0
        for line in file:
            line = line.replace("\n", "")
            list_values1.extend(map(float, line.split()[::3]))
            list_values2.extend(map(float, line.split()[1::3]))
            list_values3.extend(map(float, line.split()[2::3]))
            list_pos_xy.extend([[int(16), int(cur_deg)], [int(12), int(cur_deg)],
                                [int(8), int(cur_deg)], [int(4), int(cur_deg)]])
            cur_deg = cur_deg + 30
    list_values1 = np.array(list_values1)
    list_values2 = np.array(list_values2)
    list_values3 = np.array(list_values3)
    list_pos_xy = np.array(list_pos_xy)
    return list_values1, list_values2, list_values3, list_pos_xy


def load_function(file="Neggrid.txt"):
    list_pos_xy = []
    list_values = []
    list_pos = []
    with open(file, 'r') as file:
        cur_x = 4
        for line in file:
            cur_y = 3
            line = line.replace("\n", "")
            for value in line.split():
                list_pos.append([int(cur_x), int(cur_y), 0])
                list_pos_xy.append([int(cur_x), int(cur_y)])
                list_values.append(float(value))
                cur_y = cur_y + 1
            cur_x = cur_x + 1

    list_pos_xy = np.array(list_pos_xy)
    list_values = np.array(list_values)
    return list_values, list_pos_xy, list_pos


def correction(alpha, basis, shift, value):
    return value * (basis + (1 - basis) * ((math.cos((alpha * math.pi / 180) + shift) + 1) / 2))


def correction_rh(x, y, middle, basis, alpha, value):
    basis = abs(basis * (1 - 0.5 * (1 / middle) * abs(x - middle)))
    basis = basis * (1 / 9) * abs(9 - y)
    basis = 1 - basis
    return float(correction(alpha, basis, 0, value))


def function_shift(fileName):
    list_values, list_pos_xy, list_pos = load_function(fileName)
    aproximate_neg = interpolate.Rbf(*list_pos_xy.T, list_values, function='inverse')

    def shell_function(X_shift, Y_shift, alpha):
        if 3 <= X_shift <= 6:
            return correction_rh(X_shift, Y_shift, 4.5, 0.56, alpha, aproximate_neg(X_shift, Y_shift))
        elif 6 < X_shift <= 11:
            return correction_rh(X_shift, Y_shift, 9, 0.78, alpha, aproximate_neg(X_shift, Y_shift))
        elif 11 < X_shift:
            return correction_rh(X_shift, Y_shift, 14, 0.95, alpha, aproximate_neg(X_shift, Y_shift))

    return shell_function


def convert_cos_to_degree(cos, dist):
    if dist == "R":
        angle = 2 * 90 * cos / math.pi
    else:
        angle = 360 - 2 * 90 * cos / math.pi
    return angle


def function_shift_dipoles():
    list_values1, list_values2, list_values3, list_pos_xy = load_dipole_function("Dipgrid.txt")
    approximate_dip1 = interpolate.Rbf(*list_pos_xy.T, list_values1, function='cubic')
    approximate_dip2 = interpolate.Rbf(*list_pos_xy.T, list_values2, function='cubic')
    approximate_dip3 = interpolate.Rbf(*list_pos_xy.T, list_values3, function='cubic')

    def shell_function(X_shift, Y_shift, cos_alpha, dist="R"):
        def regression_func(x, a, b, c):
            return a * (b / x) + c

        z = []
        alpha = convert_cos_to_degree(cos_alpha, dist)
        z.append(approximate_dip1(X_shift, alpha))
        z.append(approximate_dip2(X_shift, alpha))
        z.append(approximate_dip3(X_shift, alpha))
        opt_quadratic, conv_quadratic = curve_fit(regression_func, [3.5, 4.5, 5.5], z)
        result = 0
        try:
            result = regression_func(Y_shift, *opt_quadratic)
        except ValueError:
            if Y_shift > 10:
                print(" Amino acid too far from chromophore")
                result = 0
            else:
                print("Unexpected error")
        return result

    return shell_function
