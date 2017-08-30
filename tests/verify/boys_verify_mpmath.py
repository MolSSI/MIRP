#!/usr/bin/env python3

import argparse
import sys
from mpmath import mp

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

def read_test_file(filepath):
    mtv = []
    with open(filepath, 'r') as f:
        lines = [ l for l in f.readlines() if not l.startswith('#') ] 
        bits = int(lines[0])
        for l in lines[1:]:
            m,t,v = l.strip().split(' ')
            mtv.append((m,t,v))

    return (bits, mtv)


parser = argparse.ArgumentParser()
parser.add_argument("--filename",     type=str,            required=True,  help="File to test")
parser.add_argument("--ndigits",      type=int, default=0, required=False, help="Number of decimal digits to test")
args = parser.parse_args()

test_digits, test_vals = read_test_file(args.filename)

ndigits = args.ndigits
if ndigits == 0:
    ndigits = test_digits-1

print("Verifying {} up to {} decimal digits".format(args.filename, ndigits))

if ndigits >= test_digits:
    print("!-----------------------------------------------------------------------------")
    print("! WARNING: number of digits you want to test >= number stored in the test file")
    print("!-----------------------------------------------------------------------------")

nfailed = 0
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
            nfailed += 1

            print("Bad value: {} {}".format(m, t))
            print("  In test file: {}".format(str(vtest)))
            print("        MPMath: {}".format(str(v)))
            print(" Relative diff: {}".format(str(reldiff)))
            print()

passed = ntest - nfailed
percent = 100.0*float(passed)/float(ntest)
print("{} / {} failed ({:.2f}% passed)".format(nfailed, len(test_vals), percent))
if nfailed > 0:
    quit(1)

