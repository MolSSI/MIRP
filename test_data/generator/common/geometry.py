#!/usr/bin/env python3

#
# A quick script to read XYZ files
# Also converts the coordinates to bohr
#

from mpmath import mp

def read_xyz(xyzfile_path):
    # NIST CODATA 2014 value for a bohr is used below
    ang_to_bohr = 1.0/0.52917721067

    geo = []

    with open(xyzfile_path, 'r') as f:
        # skips the first two lines
        tmp = [l.strip() for l in f.readlines()[2:]]

    # splits each line into tuples (and ignores blank lines)
    tmp = [l.split() for l in tmp if l]

    # convert to bohr
    for g in tmp:
        geo.append((g[0],
                    mp.mpf(g[1]) * ang_to_bohr,
                    mp.mpf(g[2]) * ang_to_bohr,
                    mp.mpf(g[3]) * ang_to_bohr))

    return geo
