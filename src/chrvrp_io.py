"""Instance parsing for expt1_fast.py (get_words / read_dat_file)."""

import numpy as np


def get_words(line):
    """Return a list of the tokens in line."""
    line = line.replace("\t", " ")
    line = line.replace("\v", " ")
    line = line.replace("\r", " ")
    line = line.replace("\n", " ")
    while line.count("  "):
        line = line.replace("  ", " ")
    line = line.strip()
    return [word + " " for word in line.split(" ")]


def read_dat_file(filename):
    """Parse a .tsp CHR-VRP instance file and return (s, v, d, dist).

    dist is the full N x N distance matrix, positions ordered to match the
    concatenated node ordering (suppliers, then demand, then vehicles) used
    throughout the solver. s/v/d are the supplier/vehicle/demand counts.
    """
    loc_Vx, loc_Vy = [], []
    loc_Sx, loc_Sy = [], []
    loc_Dx, loc_Dy = [], []
    with open(filename) as f:
        for line in f:
            words = get_words(line)
            start = list(words[0])[0]
            if start == 'V':
                loc_Vx.append(float(words[1]))
                loc_Vy.append(float(words[2]))
            elif start == 'S':
                loc_Sx.append(float(words[1]))
                loc_Sy.append(float(words[2]))
            elif start == 'D':
                loc_Dx.append(float(words[1]))
                loc_Dy.append(float(words[2]))

    s = len(loc_Sx)
    v = len(loc_Vx)
    d = len(loc_Dx)

    locx = loc_Sx + loc_Dx + loc_Vx
    locy = loc_Sy + loc_Dy + loc_Vy
    n = len(locx)

    dist = []
    for i in range(n):
        dist.append([])
        x1 = locx[i]
        y1 = locy[i]
        for j in range(n):
            d1 = (x1 - locx[j]) * (x1 - locx[j])
            d2 = (y1 - locy[j]) * (y1 - locy[j])
            dist[i].append(np.sqrt(d1 + d2) / 1000)

    return s, v, d, dist
