#!/usr/bin/env python

################################################################################
# Analyzing the pi-orbital axis vector and sigma-pi hybridization based on     #
# geometry written in .gui file (CRYSTAL fort.34 format).                      #
# Output is written in the atomic data format of latt_plot.py                  #
# Used only for Carbon atoms. For theoretical background, please refer to:     #
# DOI: 10.1021/jp010023f                                                       #
# ----------------------------                                                 #
# by Spica. Vir., ICL, Aug. 04 - 21                                            #
# spica.h.zhou@gmail.com                                                       #
################################################################################

import numpy as np
import math
import sys
import copy


class Bond:
    """
    Class bond, used only for connectivity analysis. No data is written in the 
    data output.
    """

    def __init__(self, data=None, 
                 bg_label=0, bg_fx=0., bg_fy=0.,
                 ed_label=0, ed_fx=1., ed_fy=1.):
        self.data = data
        self.bg_label = bg_label
        self.ed_label = ed_label
        self.bg_fx = bg_fx
        self.bg_fy = bg_fy
        self.ed_fx = ed_fx
        self.ed_fy = ed_fy

    def write_data(self, outfile):
        outfile.write('%-5d%5d%12.8f\n' % (self.bg_label,
                                           self.ed_label,
                                           self.data))


class Atom:
    """
    Class atom, stands for a single atom, with information of label and 
    coordinate. 
    """

    def __init__(self, label=0, atomic_num=0, data=[],
                 x=0., y=0., z=0.,
                 fracx=0., fracy=0., bonded_label=[]):
        self.label = label
        self.atomic_num = atomic_num
        self.data = data
        self.x = x
        self.y = y
        self.z = z
        self.fracx = fracx
        self.fracy = fracy
        self.bonded_label = bonded_label

    def write_data(self, outfile, seq=1):
        if seq == 1:
            if self.data[0] != None:
                outfile.write('%-8d%12.8f\n' % (self.label, self.data[0]))
            else:
                return
        else: 
            if self.data[seq - 1] != None:
                outfile.write('%-8d%12.8f\n' % (self.label, self.data[seq - 1]))
            else:
                return


class Cell:
    """
    Define the lattice for calculation, i.e., the lattice defined in fort.34
    file. Plotted only if an expanding matrix different from default value is 
    defined. 

    vect_a/b includes the 2D lattice vectors of the cell
    """

    def __init__(self, vect_a=[], vect_b=[]):
        self.vect_a = vect_a
        self.vect_b = vect_b


class Pyramid:
    """
    constructing Pyramid cluster and get translation matrix between cartesian 
    coordinates and the internal coordinates defined by the central atom and its 
    3 neighbors. 

    POAV and orbital hybridization are calculated using the normalized internal 
    coordinate. Theoretical basis: 10.1021/jp010023f
    """

    def __init__(self, cell, central_atom_lb=0, atom_list=[], allowed_dist=2.0):
        self.c_atom = atom_list[central_atom_lb - 1]

        near_atom = []
        for i in atom_list[central_atom_lb - 1].bonded_label:
            near_atom.append(atom_list[i - 1])
        
        near_atom_coord = []
        for i in near_atom:
            for a in range(-1, 2):
                for b in range(-1, 2):
                    xcoord = i.x + a * cell.vect_a[0] + b * cell.vect_b[0]
                    ycoord = i.y + a * cell.vect_a[1] + b * cell.vect_b[1]
                    zcoord = i.z
                    dist = ((self.c_atom.x - xcoord) ** 2 +
                            (self.c_atom.y - ycoord) ** 2 +
                            (self.c_atom.z - zcoord) ** 2)
                    if dist <= allowed_dist ** 2:
                        break

                if dist <= allowed_dist ** 2:
                        break

            near_atom_coord.append([xcoord, ycoord, zcoord])

        self.n_atom = near_atom_coord

    def poav_anal(self):
        if len(self.n_atom) != 3:
            print('Warning: Atom ' + str(self.c_atom.label) 
                + ' the number of sigma bonds should be 3. No output data here.')
            poav = None
            sig_pi_hybr = None
            return poav, sig_pi_hybr

        central = [self.c_atom.x, self.c_atom.y, self.c_atom.z]
        central = np.array(central, dtype = float)
        near1 = np.array(self.n_atom[0], dtype = float)
        near2 = np.array(self.n_atom[1], dtype = float)
        near3 = np.array(self.n_atom[2], dtype = float)

        vect1 = near1 - central
        vect1 = vect1 / np.linalg.norm(vect1)
        vect2 = near2 - central
        vect2 = vect2 / np.linalg.norm(vect2)
        vect3 = near3 - central
        vect3 = vect3 / np.linalg.norm(vect3)

        new_x = -vect1
        new_z = np.cross(vect1, vect2)
        new_z = new_z / np.linalg.norm(new_z)
        new_y = np.cross(new_z, new_x)
        new_y = new_y / np.linalg.norm(new_y)

        internal = np.matrix([new_x, new_y, new_z])
        internal_rev = np.linalg.inv(internal)

        vect1_int = np.array(vect1 * internal_rev)[0]
        vect2_int = np.array(vect2 * internal_rev)[0]
        vect3_int = np.array(vect3 * internal_rev)[0]
        # print(self.c_atom.label)
        # print([central[0], central[1], central[2]])
        # print([near1[0], near1[1], near1[2]]) 
        # print([near2[0], near2[1], near2[2]])
        # print([near3[0], near3[1], near3[2]])

        base = (vect2_int[1] * vect3_int[2]) ** 2 + \
            (vect3_int[2] * (vect2_int[0] - vect1_int[0])) ** 2 + \
            (vect3_int[1] * (vect2_int[0] - vect1_int[0]) - \
                vect2_int[1] * (vect3_int[0] - vect1_int[0])) ** 2 
        cos_sigmapi = (vect1_int[0] * vect2_int[1] * vect3_int[2]) / (base ** 0.5)
        sigmapi = math.acos(cos_sigmapi) / np.pi * 180
        poav = round(sigmapi - 90, 8)

        sig_pi_hybr = round(2 * cos_sigmapi ** 2 / (1 - 3 * cos_sigmapi ** 2), 8)

        return poav, sig_pi_hybr


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
            latt_cell = Cell(vect_a=[float(vect_a[0]), float(vect_a[1])],
                             vect_b=[float(vect_b[0]), float(vect_b[1])])
            latvec_matrix = np.matrix([latt_cell.vect_a, latt_cell.vect_b],
                                      dtype=float)
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
                frac = np.array(
                    [float(coord[1]), float(coord[2])]) * latvec_rev
                if frac[0, 0] < 0.:
                    frac[0, 0] += 1

                if frac[0, 1] < 0.:
                    frac[0, 1] += 1

                xycoord = frac * latvec_matrix
                an_atom = Atom(label=int(countlabel),
                               atomic_num=int(coord[0]),
                               data=[],
                               x=float(xycoord[0, 0]),
                               y=float(xycoord[0, 1]), z=float(coord[3]),
                               fracx=frac[0, 0], fracy=frac[0, 1],
                               bonded_label=[])
                atoms.append(an_atom)
            else:
                natom = int(gui_line[countline])
                countline += 1
                continue

        countline += 1

    return atoms, natom, latt_cell


def connectivity(atoms, natom, latt_cell, bond_max=2.0):
    """
    The same algorithm as in latt_plot.py is adopted.
    """
    atoms_new, natom_new = cellexpand(latt_cell, atoms=atoms,
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
                dist = ((i.x - atoms_new[j].x) ** 2 + (i.y - atoms_new[j].y)
                        ** 2 + (i.z - atoms_new[j].z) ** 2) ** 0.5
                if dist <= bond_max:
                    if rep == 0:
                        bonded_label.append(label_judge)

                    a_bond = Bond(bg_label=i.label,
                                  bg_fx=i.fracx, bg_fy=i.fracy,
                                  ed_label=atoms_new[j].label,
                                  ed_fx=atoms_new[j].fracx,
                                  ed_fy=atoms_new[j].fracy)
                    bonds.append(a_bond)
                    i.bonded_label.append(num)
                    atoms[num - 1].bonded_label.append(i.label)

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
                        xycoord = np.dot(
                            np.matrix([fracx, fracy]), latvec_matrix)
                        an_atom = Atom(label=i.label,
                                       atomic_num=i.atomic_num,
                                       x=xycoord[0, 0], y=xycoord[0, 1],
                                       z=i.z, data=i.data,
                                       fracx=fracx, fracy=fracy,
                                       bonded_label=[])
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
                        a_bond = Bond(bg_label=i.bg_label,
                                      bg_fx=bg_fx, bg_fy=bg_fy,
                                      ed_label=i.ed_label,
                                      ed_fx=ed_fx, ed_fy=ed_fy)
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

for i in atoms:
    poav_calculator = Pyramid(cell = latt_cell, 
                              central_atom_lb = i.label, 
                              atom_list = atoms)
    poav, sig_pi_hybr = poav_calculator.poav_anal()
    i.data = [poav, sig_pi_hybr]

wtatom1 = open("atom_poav.dat", "w")
wtatom1.write('%-8s%12s\n' % ('LABEL', 'POAV (deg)'))
for i in atoms:
    i.write_data(outfile = wtatom1, seq = 1)

wtatom1.close()

wtatom2 = open("atom_hybridisation.dat", "w")
wtatom2.write('%-8s%12s\n' % ('LABEL', 'S-P HYBRID'))
for i in atoms:
    i.write_data(outfile = wtatom2, seq = 2)

wtatom2.close()
