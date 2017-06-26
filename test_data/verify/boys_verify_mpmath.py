#!/usr/bin/env python3

import argparse
import sys

from mpmath import mp

import mirppy
import mirppy.testing_utilities

def calculate_boys(m, t):
    zero_mp = mp.mpf(0)
    one_mp = mp.mpf(1)
    two_mp = mp.mpf(2)
    m_mp = mp.mpf(m)
    t_mp = mp.mpf(t)

    if t_mp == zero_mp:
        F_mp = one_mp/(two_mp * m_mp + one_mp)
    else:
        M = m_mp + (one_mp / two_mp)
        F_mp = mp.gammainc(M, zero_mp, t_mp) * one_mp/(two_mp * mp.power(t_mp, M))
    return F_mp


parser = argparse.ArgumentParser()
parser.add_argument("--filename",     type=str,            required=True,  help="File to test")
parser.add_argument("--ndigits",      type=int, default=0, required=False, help="Number of decimal digits to test")
args = parser.parse_args()

test_digits, test_vals = mirppy.testing_utilities.read_test_file(args.filename)

ndigits = args.ndigits
if ndigits == 0:
    ndigits = test_digits-1

bad_vals = []
passed = 0
ntest = len(test_vals)

with mp.workdps(test_digits+8):
    test_eps = mp.power(10, -ndigits)

    for m,t,vtest in test_vals:
        m = int(m)
        v = calculate_boys(m, mp.mpf(t))

        v = mp.mpf(v)
        vtest = mp.mpf(vtest)

        same = mp.almosteq(v, vtest, abs_eps=0, rel_eps=test_eps)

        if not same:
            reldiff = (mp.fabs(v-vtest))/max(mp.fabs(v), mp.fabs(vtest))
            bad_vals.append((m, t, str(vtest), str(v), reldiff))
        else:
            passed += 1

    for m,t,vtest,v,reldiff in bad_vals:
        print("Bad value: {} {}".format(m, t))
        print("  In test file: {}".format(str(vtest)))
        print("        MPMath: {}".format(str(v)))
        print(" Relative diff: {}".format(str(reldiff)))
        print()

percent = 100.0*float(passed)/float(ntest)
print("{} / {} passed ({:.2f}%)".format(passed, len(test_vals), percent))
if len(bad_vals) > 0:
    quit(1)

