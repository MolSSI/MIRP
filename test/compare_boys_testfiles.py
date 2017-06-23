import argparse
from mpmath import mp

import mirppy
import mirppy.testing_utilities

parser = argparse.ArgumentParser()
parser.add_argument("--filename", type=str, required=True, help="Output file name")
args = parser.parse_args()

prec, mtvdata = mirppy.testing_utilities.read_test_file(args.filename)


failed_interval = []

with mp.workprec(prec):
    eps_compare = mp.power(2, -prec+3)
    for m,t,v in mtvdata: 
        myval = mirppy.mirp_boys_interval(m, t, prec*2)[m]
        v1 = mp.mpf(v)
        v2 = mp.mpf(myval)
        if not mp.almosteq(v1, v2, rel_eps=eps_compare):
            failed_interval.append((m, t, v1, v2))

print("{} / {} Failed Testing".format(len(failed_interval), len(mtvdata)))
with mp.workprec(prec):
    for i in failed_interval:
        print("{} {}".format(i[0], i[1]))
        print(i[2])
        print(i[3])
        print((i[2]-i[3]))
        print((i[2]-i[3])/i[2])
        print()

if len(failed_interval) > 0:
    quit(1)

