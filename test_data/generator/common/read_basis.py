#!/usr/bin/env python3

#
# A quick script to read an XYZ file
# Also converts the coordinates to bohr
#

from mpmath import mp

def read_basis(basfile_path):
    # Read basis set
    with open(basfile_path, 'r') as f:
        flines = [ l.strip() for l in f.readlines() ]

    # Strip out comments and blank lines
    flines = [ l for l in flines if l and not l.startswith('!') ]

    # atom -> basis mapping
    atombas = { }

    i = 1  # skip initial **** separator
    while i < len(flines):
        atom_symbol = flines[i].split()[0]
        atom = []

        i += 1
        while flines[i] != '****':
            lsplt = flines[i].split()
            stype = lsplt[0].upper()
            nshell = int(lsplt[1])

            alpha = []
            coeff = []
            i += 1
            for j in range(0, nshell):
                lsplt = flines[i].split()
                alpha.append(mp.mpf(lsplt[0]))
                coeff.append([ mp.mpf(c) for c in lsplt[1:] ])
                i += 1


            nprim = len(alpha)

            # Transposes the coefficient matrix
            coeff = list(map(list, zip(*coeff)))

            # Handle sp, spd, etc
            if len(stype) > 1:
                for j,s in enumerate(stype): 
                    shell = { 'type' : s,
                              'alpha' : alpha,
                              'coeff' : coeff[j],
                              'nprim' : nprim,
                              'ngeneral' : 1
                            }

                    atom.append(shell)
                    
            else:
                shell = { 'type' : stype,
                          'alpha' : alpha,
                          'coeff' : coeff[0],
                          'nprim' : nprim,
                          'ngeneral' : len(coeff[0]) // nprim
                        }

                atom.append(shell)

        i += 1

        # add atom to basis set
        atombas[atom_symbol] = atom


    return atombas

#    # Apply basis set to molecule
#    mol = []
#
#    for a in geo:
#      mol.append({
#                   'Sym' : a[0],
#                   'XYZ' : [ float(a[1])*ang_to_bohr, float(a[2])*ang_to_bohr, float(a[3])*ang_to_bohr ],
#                   'BAS' : atombas[a[0]]
#                 })
#
#
#    print("-----------------------------------")
#
#    # Create file
#    with open(args.out, 'w') as f:
#      f.write(str(len(mol)) + "\n")
#      for a in mol:
#        f.write("{} {} {} {}\n".format(a['Sym'], len(a['BAS']['Shells']), a['BAS']['NPRIM1'], a['BAS']['NPRIM2']))
#        f.write("{} {} {}\n".format(a['XYZ'][0], a['XYZ'][1], a['XYZ'][2]))
#
#        for b in a['BAS']['Shells']:
#          nprim = len(b['Alpha'])
#          ngen = len(b['Coef'][0])
#          f.write('{} {} {}\n'.format(b['Type'], nprim, ngen))
#
#          for i in range(0, nprim):
#            f.write('{}'.format(b['Alpha'][i]))
#
#            for j in range(0, ngen):
#              f.write('    {}'.format(b['Coef'][i][j]))
#
#            f.write("\n")
#
