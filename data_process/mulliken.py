#!/usr/bin/env python

################################################################################
# To process mulliken charge output of CRYSTAL 17 (PPAN)                       #
# The total atomic charge of a specific atom is written in the input form of   #
# latt_plot.py                                                                 #
# Note: unsuitable for atoms with only one atomic orbital                      #
# ----------------------------                                                 #
# by Spica. Vir., ICL, Aug. 03 - 21                                            #
# spica.h.zhou@gmail.com                                                       #
################################################################################

import sys
import math


class Atom:
    """
    storing data of atom labels and atomic charge
    """

    def __init__(self, label=0, charge=[], orbit_chg=[]):
        self.label = label
        self.charge = charge
        self.orbit_chg = orbit_chg

    def writedata(self, output, spec_chg=False):
        if not spec_chg:
            for i in self.charge:
                output.write('%-8d%12.4f\n' % (self.label, i))

        else:
            for i in self.orbit_chg:
                output.write('%-8d%12.4f\n' % (self.label, i))


def read_data(ppan_line, atoms, natom, specific_orbit, bg_line=8):
    """
    Read the total atomic charge. Note the sequence of atoms should be kept
    exactly in the same way as output geometry file (.f34 or .gui), so the raw
    output without any modification is preferred.
    """
    # skip notification lines / spin-positive lines
    countline = bg_line
    countatom = 0
    nshell = int(ppan_line[bg_line].strip().split()[1])

    while countatom < natom:
        nshell = int(ppan_line[countline].strip().split()[1])
        countline += 2
        atomic_chg = float(ppan_line[countline].strip().split()[0])
        atoms[countatom].label = countatom + 1
        atoms[countatom].charge.append(atomic_chg)
        locat_norbit = (nshell + 1) // 8 + 1
        countline += locat_norbit
        norbit = int(ppan_line[countline].strip().split()[0])
        norbitlines = math.ceil(norbit / 8)
        orbit_chg = []
        if specific_orbit:
            for i in range(norbitlines):
                countline += 1
                orbit_chg += ppan_line[countline].strip().split()

            ochg = 0
            for j in specific_orbit:
                ochg += float(orbit_chg[j - 1])

            atoms[countatom].orbit_chg.append(ochg)
        else:
            countline += norbitlines

        countline += 1
        countatom += 1

    return atoms, countline


# main I/O function
ppan_file = sys.argv[1]
ppan = open(ppan_file, "r")
ppan_line = ppan.readlines()
ppan.close()
natom = int(ppan_line[7].strip().split()[1])
nspin = int(ppan_line[7].strip().split()[0])
######################### modify accords to needs ##############################
# specific orbit for atomic charge or spin population is kept the same. 
# e.g., C atom pz orbital, 6-21G* Basis set. (Also, refer to output for numbers)
specific_orbit = [5, 9]
# specific_orbit = []
################################################################################
atoms = []
for i in range(natom):
    an_atom = Atom(label = 0, charge = [], orbit_chg = [])
    atoms.append(an_atom)

atoms, countline = read_data(ppan_line, atoms, natom, specific_orbit)

if nspin == 2:
    atoms_spin = []
    for i in range(natom):
        an_atom_spin = Atom(label = 0, charge = [], orbit_chg = [])
        atoms_spin.append(an_atom_spin)

    atoms_spin, countline = read_data(ppan_line, atoms_spin, natom, 
                                      specific_orbit, bg_line = countline)

output = open("mulliken_total.dat", "w")
output.write('%-8s%12s\n' % ('LABEL', 'Tot_chg'))
for i in atoms:
    i.writedata(output = output)

output.close()

if specific_orbit:
    output_specific = open("mulliken_specific.dat", "w")
    output_specific.write('%-8s%12s\n' % ('LABEL', 'Spec_chg'))
    for i in atoms:
        i.writedata(output = output_specific, spec_chg = True)

    output_specific.close()

if nspin == 2:
    output_spin = open("mulliken_spin.dat", "w")
    output_spin.write('%-8s%12s\n' % ('LABEL', 'Spin_tot'))
    for i in atoms_spin:
        i.writedata(output = output_spin)

    output_spin.close()

    if specific_orbit:
        output_spin_specific = open("mulliken_spin_specific.dat", "w")
        output_spin_specific.write('%-8s%12s\n' % ('LABEL', 'Spin_spec'))
        for i in atoms_spin:
            i.writedata(output = output_spin_specific, spec_chg = True)

        output_spin_specific.close()
