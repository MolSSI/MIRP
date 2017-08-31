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
    flines = [l for l in flines if l and not l.startswith('#') and not l.startswith('BASIS')]

    # atom -> basis mapping
    atom_basis = {}

    i = 0
    while flines[i] != "END":
        lsplt = flines[i].split()
        stype = lsplt[1].lower()
        atom_symbol = lsplt[0]

        i += 1

        alpha = []
        coeff = []
        while flines[i][0].isdigit():
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

                if not atom_symbol in atom_basis:
                    atom_basis[atom_symbol] = [shell]
                else:
                    atom_basis[atom_symbol].append(shell)

        else:
            shell = {'am': amchar_to_int(stype),
                     'alpha': alpha,
                     'coeff': coeff[0],
                     'nprim': nprim,
                     'ngeneral': len(coeff[0]) // nprim
                     }

            if not atom_symbol in atom_basis:
                atom_basis[atom_symbol] = [shell]
            else:
                atom_basis[atom_symbol].append(shell)

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
