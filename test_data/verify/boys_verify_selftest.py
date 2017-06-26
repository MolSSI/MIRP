#!/usr/bin/env python3

import argparse
import sys

from mpmath import mp

import mirppy
import mirppy.testing_utilities

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
        v1 = mirppy.mirp_boys_createtest(m, t, test_digits)[m]
        v2 = mirppy.mirp_boys_interval(m, t, mp.prec)[m]
        v3 = mirppy.mirp_boys_mp(m, t, mp.prec)[m]

        v1 = mp.mpf(v1)
        v2 = mp.mpf(v2)
        v3 = mp.mpf(v3)
        vtest = mp.mpf(vtest)

        same1 = mp.almosteq(v1, vtest, abs_eps=0, rel_eps=test_eps)
        same2 = mp.almosteq(v2, vtest, abs_eps=0, rel_eps=test_eps)
        same3 = mp.almosteq(v3, vtest, abs_eps=0, rel_eps=test_eps)

        if not (same1 and same2 and same3):
            bad_vals.append((m, t, str(vtest), str(v1), str(v2), str(v3)))
        else:
            passed += 1

    for m,t,v1,v2,v3,vtest in bad_vals:
        print("Bad value: {} {}".format(m, t))
        print("  In test file: {}".format(str(vtest)))
        print("    CreateTest: {}".format(str(v1)))
        print("      Interval: {}".format(str(v2)))
        print("          MPFR: {}".format(str(v3)))
        print()

percent = 100.0*float(passed)/float(ntest)
print("{} / {} passed ({:.2f}%)".format(passed, len(test_vals), percent))
if len(bad_vals) > 0:
    quit(1)

