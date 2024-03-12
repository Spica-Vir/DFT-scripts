#!/usr/bin/env python

import numpy as np
import math
import sys
import copy

class Bond:
    """
        Class bond, plotted as colored/grey lines.
    """

    def __init__(self, data=None, x=[0., 1.], y=[0., 1.], 
                 bg_label=0, bg_fx=0., bg_fy=0., 
                 ed_label=0, ed_fx=1., ed_fy=1., ):
        self.data = data
        self.x = x
        self.y = y
        self.bg_label = bg_label
        self.ed_label = ed_label
        self.bg_fx = bg_fx
        self.bg_fy = bg_fy
        self.ed_fx = ed_fx
        self.ed_fy = ed_fy

    def write_data(self, outfile):
        outfile.write('%-5d%5d%12.8f\n' %(self.bg_label, 
                                          self.ed_label, 
                                          self.data))


class Atom:
    """
        Class atom, plotted as colored circles. 
    """

    def __init__(self, label=0, atomic_num=0, data=None, 
                 x=0., y=0., z=0., 
                 fracx=0., fracy=0.):
        self.label = label
        self.atomic_num = atomic_num
        self.data = data
        self.x = x
        self.y = y
        self.z = z
        self.fracx = fracx
        self.fracy = fracy
    
    def write_data(self, outfile):
        outfile.write('%-8d%12.8f\n' %(self.label, self.data))


class Cell:
    """
        Define the lattice for calculation, i.e., the lattice defined in fort.34
        file. Plotted only if an expanding matrix different from default value 
        is defined. 

        vect_a/b includes the 2D lattice vectors of the cell
    """

    def __init__(self, vect_a=[], vect_b=[]):
        self.vect_a = vect_a
        self.vect_b = vect_b

def read_geom(gui_line):
    """
        read the geometry stored in .gui file. 
        atoms is the list of Atom class, with coordinates and labels of atoms
        natom is the number of atoms
        latt_cell is the Cell class, with list of cell edges. 
    """
    countline = 0
    nsymm = 0
    natom = 0
    atoms = []
    countlabel = 0
    while countline < len(gui_line) - 1:
        if countline == 1:
            vect_a = gui_line[countline].strip().split()
            vect_b = gui_line[countline + 1].strip().split()
            latt_cell = Cell(vect_a = [float(vect_a[0]), float(vect_a[1])],
                             vect_b = [float(vect_b[0]), float(vect_b[1])])
            latvec_matrix = np.matrix([latt_cell.vect_a, latt_cell.vect_b], 
                                      dtype = float)
            latvec_rev = np.linalg.inv(latvec_matrix)
            countline += 3
            continue

        if countline == 4:
            nsymm = int(gui_line[countline])
            countline += int(nsymm * 4) + 1
            continue

        if nsymm != 0:
            if natom != 0:
                countlabel += 1
                coord = gui_line[countline].strip().split()
                frac = np.array([float(coord[1]), float(coord[2])]) * latvec_rev
                if frac[0, 0] < 0.:
                    frac[0, 0] += 1

                if frac[0, 1] < 0.:
                    frac[0, 1] += 1

                xycoord = frac * latvec_matrix
                an_atom = Atom(label = int(countlabel), 
                               atomic_num = int(coord[0]), 
                               data = float(coord[3]), 
                               x = float(xycoord[0, 0]), 
                               y = float(xycoord[0, 1]), z = float(coord[3]), 
                               fracx = frac[0, 0], fracy = frac[0, 1])
                atoms.append(an_atom)
            else:
                natom = int(gui_line[countline])
                countline += 1
                continue

        countline += 1

    return atoms, natom, latt_cell

def connectivity(atoms, natom, latt_cell, bond_max=2.0):
    atoms_new, natom_new = cellexpand(latt_cell, atoms = atoms,
                                      expand=[[-1, 1], [0, 2]])
    bonds = []
    bonded_label = []
    for i in atoms:
        for j in range(natom_new):
            num = atoms_new[j].label
            label_judge = [min([i.label, num]), max([i.label, num])]

            if math.isclose(atoms[num - 1].x, atoms_new[j].x, rel_tol=1e-2) & \
               math.isclose(atoms[num - 1].y, atoms_new[j].y, rel_tol=1e-2): 
                rep = 0
            else:
                rep = 1

            if ((rep != 0) & (i.label != num)) | \
            ((rep == 0) & (i.label != num) & (not label_judge in bonded_label)):
                dist = ((i.x - atoms_new[j].x) ** 2 + (i.y - atoms_new[j].y)\
                 ** 2 + (i.z - atoms_new[j].z) ** 2) ** 0.5
                if dist <= bond_max:
                    if rep == 0:
                        bonded_label.append(label_judge)

                    a_bond = Bond(data = dist, 
                                  x = [i.x, atoms_new[j].x], 
                                  y = [i.y, atoms_new[j].y], 
                                  bg_label = i.label, 
                                  bg_fx = i.fracx, bg_fy = i.fracy, 
                                  ed_label = atoms_new[j].label,
                                  ed_fx = atoms_new[j].fracx, 
                                  ed_fy = atoms_new[j].fracy)
                    bonds.append(a_bond)

    return bonds

def cellexpand(latt_cell, atoms=[], bonds=[], expand=[[0, 1], [0, 1]]):
    origin_xmi = math.floor(expand[0][0])
    origin_xmx = math.ceil(expand[0][1])
    origin_ymi = math.floor(expand[1][0])
    origin_ymx = math.ceil(expand[1][1])

    latvec_matrix = np.matrix([latt_cell.vect_a, latt_cell.vect_b])
    natom = len(atoms)

    if atoms and not bonds:
        if expand == [[0, 1], [0, 1]]:
            return atoms, natom

        atoms_new = []
        for x in range(origin_xmi, origin_xmx):
            for y in range(origin_ymi, origin_ymx):
                for i in atoms: 
                    fracx = i.fracx + x
                    fracy = i.fracy + y
                    if (fracx > expand[0][0]) & (fracx < expand[0][1]) & \
                    (fracy > expand[1][0]) & (fracy < expand[1][1]):
                        xycoord = np.dot(np.matrix([fracx, fracy]), latvec_matrix)
                        an_atom = Atom(label = i.label, 
                                       atomic_num = i.atomic_num, 
                                       x = xycoord[0, 0], y = xycoord[0, 1],  
                                       z = i.z, data = i.data, 
                                       fracx = fracx, fracy = fracy)
                        atoms_new.append(an_atom)

        natom_new = len(atoms_new)
        return atoms_new, natom_new

    if bonds and not atoms:
        if expand == [[0, 1], [0, 1]]:
            return bonds

        bonds_new = []
        for x in range(origin_xmi, origin_xmx):
            for y in range(origin_ymi, origin_ymx):
                for i in bonds:
                    bg_fx = i.bg_fx + x
                    bg_fy = i.bg_fy + y
                    ed_fx = i.ed_fx + x
                    ed_fy = i.ed_fy + y
                    judge_bg = (bg_fx > expand[0][0]) & (bg_fx < expand[0][1]) \
                             & (bg_fy > expand[1][0]) & (bg_fy < expand[1][1])
                    judge_ed = (ed_fx > expand[0][0]) & (ed_fx < expand[0][1]) \
                             & (ed_fy > expand[1][0]) & (ed_fy < expand[1][1])
                    if judge_bg | judge_ed:
                        xy = np.dot(np.matrix([[bg_fx, bg_fy], [ed_fx, ed_fy]]), 
                                              latvec_matrix)
                        a_bond = Bond(data = i.data, color = i.color, 
                                      x = [xy[0, 0], xy[1, 0]], 
                                      y = [xy[0, 1], xy[1, 1]], 
                                      bg_label = i.bg_label, 
                                      bg_fx = bg_fx, bg_fy = bg_fy, 
                                      ed_label = i.ed_label, 
                                      ed_fx = ed_fx, ed_fy = ed_fy)
                        bonds_new.append(a_bond)

        return bonds_new

    if bonds and atoms:
        return atoms, bonds
# main function
gui_file = sys.argv[1]
gui = open(gui_file, "r")
gui_line = gui.readlines()
gui.close()

atoms, natom, latt_cell = read_geom(gui_line)
bonds = connectivity(atoms, natom, latt_cell)

wtatom = open("atom_height.dat", "w")
wtatom.write('%-8s%12s\n' %('LABEL', 'HEIGHT (A)'))
for i in atoms:
    i.write_data(outfile = wtatom)

wtatom.close()

wtbond = open("bond_length.dat", "w")
wtbond.write('%-5s%5s%12s\n' %('BG', 'ED', 'LENGTH (A)'))
for i in bonds:
    i.write_data(outfile = wtbond)

wtbond.close()
