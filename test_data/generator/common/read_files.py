#!/usr/bin/env python3

#
# A quick script to read an XYZ file
# Also converts the coordinates to bohr
#

from mpmath import mp


def amchar_to_int(amchar):
    lookup = 'spdfghijklmnoqrtuvwxyzabce'
    return lookup.index(amchar)


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


def read_basis(basisfile_path):
    # Read basis set
    with open(basisfile_path, 'r') as f:
        flines = [l.strip() for l in f.readlines()]

    # Strip out comments and blank lines
    flines = [l for l in flines if l and not l.startswith('!')]

    # atom -> basis mapping
    atom_basis = {}

    i = 1  # skip initial **** separator
    while i < len(flines):
        atom_symbol = flines[i].split()[0]
        atom = []

        i += 1
        while flines[i] != '****':
            lsplt = flines[i].split()
            stype = lsplt[0].lower()
            nshell = int(lsplt[1])

            alpha = []
            coeff = []
            i += 1
            for j in range(0, nshell):
                lsplt = flines[i].split()
                alpha.append(mp.mpf(lsplt[0]))
                coeff.append([mp.mpf(c) for c in lsplt[1:]])
                i += 1

            nprim = len(alpha)

            # Transposes the coefficient matrix
            coeff = list(map(list, zip(*coeff)))

            # Handle sp, spd, etc
            if len(stype) > 1:
                for j, s in enumerate(stype):
                    shell = {'am': amchar_to_int(s),
                             'alpha': alpha,
                             'coeff': coeff[j],
                             'nprim': nprim,
                             'ngeneral': 1
                             }

                    atom.append(shell)

            else:
                shell = {'am': amchar_to_int(stype),
                         'alpha': alpha,
                         'coeff': coeff[0],
                         'nprim': nprim,
                         'ngeneral': len(coeff[0]) // nprim
                         }

                atom.append(shell)

        i += 1

        # add atom to basis set
        atom_basis[atom_symbol] = atom

    return atom_basis


def apply_basis(geometry, atombasis):
    shells = []

    for atom in geometry:
        b = atombasis[atom[0]]
        for bshell in b:
            shells.append(((atom[1], atom[2], atom[3]), bshell))
    return shells


def construct_basis(xyzfile_path, basisfile_path):
    geo = read_xyz(xyzfile_path)
    bas = read_basis(basisfile_path)
    return apply_basis(geo, bas)
