#!/usr/bin/env python3

from mpmath import mp
from .geometry import read_xyz 

def amchar_to_int(amchar):
    lookup = 'spdfghijklmnoqrtuvwxyzabce'
    return lookup.index(amchar)


def read_basis(basisfile_path):
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
            shells.append((tuple(atom[:4]), bshell))
    return shells


def construct_basis(xyzfile_path, basisfile_path):
    geo = read_xyz(xyzfile_path)
    bas = read_basis(basisfile_path)
    mol = apply_basis(geo, bas)
    return mol
       
 
def uncontract_basis(basis):
    shells = []
    for shell in basis:
        for a in shell[1]['alpha']:
            shells.append((shell[0], { 'am': shell[1]['am'],
                                       'alpha': [a],
                                       'coeff': [mp.mpf("1.0")],
                                       'nprim': 1,
                                       'ngeneral' : 1
                                      }))
    return shells
