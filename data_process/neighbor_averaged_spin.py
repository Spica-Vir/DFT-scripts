#!/usr/bin/env python

################################################################################
# To calculate the average mulliken spin population of the ith nearest         #
# neighbor of the central atom. Geometry input read from .gui file (CRYSTAL    #
# fort.34 format). Mulliken spin population read from script mulliken.py.      #
# Note: supercell expansion not suggested in this script, because the center   #
# of interest will be periodic. Therefore connectivity across cell boundaries  #
# will not be analyzed.                                                        #
# ----------------------------                                                 #
# by Spica. Vir., ICL, Aug. 11 - 21                                            #
# spica.h.zhou@gmail.com                                                       #
################################################################################

import numpy as np
import math
import sys
import copy


class Atom:
    """
    Class atom, stands for a specific atom in cell. 
    """

    def __init__(self, label=0, atomic_num=0, data=None,
                 x=0., y=0., z=0., bonded_atom=[]):
        self.label = label
        self.atomic_num = atomic_num
        self.data = data
        self.x = x
        self.y = y
        self.z = z
        self.bonded_atom = bonded_atom

class Cell:
    """
    Define the lattice for coordinate transformation. 

    vect_a/b includes the 2D lattice vectors of the cell
    """

    def __init__(self, vect_a=[], vect_b=[]):
        self.vect_a = vect_a
        self.vect_b = vect_b


class Neighbor_cluster:
    """
    build the cluster of atoms (the ith nearest atoms of the central atom) and 
    calculate the average Mulliken spin population of the ith neighbor. 
    """

    def __init__(self, central_atom_lb=0, neighbor_seq=0):
        self.c_atom = central_atom_lb
        self.n_seq = neighbor_seq
        

    def findneighbor(self, atom_list):
        """
        Find the labels of ith nearest neighboring atoms. Repeat count conunt is
        eliminated within the same layer. Iternations scanning back and forth 
        near central atom is also avoided. However, there still exists repeat 
        scanning over the same atom (i.e., an atom belongs to differnet layers)
        after initial layers. This problem is solved in method average_spin. 
        """
        if self.n_seq == 0:
            self.second_lb = [[self.c_atom]]
            return self

        first_bond = [self.c_atom]
        neighbors = [[self.c_atom]]
        second_bond = []
        for i in range(self.n_seq):
            moveback = max(i - 1, 0)
            for j in first_bond:
                for a in atom_list[j - 1].bonded_atom:
                    if a not in neighbors[moveback]:
                        lista = [a]
                        second_bond = second_bond + lista
                    else:
                        continue

            second_bond = list(set(second_bond))
            neighbors.append(second_bond)
            first_bond = second_bond
            second_bond = []

        self.second_lb = neighbors

        return self

    def average_spin(self, atom_list, just_this_neighbor=True):
        """
        Compute the averaged spin density of ith neighbor. If just_this_neighbor
        equals Flase, code calculates the averaged spin density from the central
        atom. Otherwise only data from this layer will be calculated. 
        """
        total_spin = 0
        total_atom = 0
        if just_this_neighbor:
            if len(self.second_lb) > 1:
                neighbor = self.second_lb[-1]
                others = self.second_lb[0:-1]
                others = [j for i in others for j in i]
                for i in neighbor:
                    if i not in others:
                        total_spin += atom_list[i - 1].data
                        total_atom += 1

            else: 
                total_spin += atom_list[self.second_lb[0][0] - 1].data
                total_atom += 1

        else:
            self.second_lb = [j for i in self.second_lb for j in i]
            self.second_lb = list(set(self.second_lb))
            for i in self.second_lb:
                total_spin += atom_list[i - 1].data
                total_atom += 1

        mean_spin = total_spin / total_atom
        # print(total_atom)

        return mean_spin


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
                an_atom = Atom(label = int(countlabel), 
                               atomic_num = int(coord[0]), data = None, 
                               x = float(xycoord[0, 0]), y = float(xycoord[0, 1]), 
                               z = float(coord[3]), bonded_atom = [])
                atoms.append(an_atom)
            else:
                natom = int(gui_line[countline])
                countline += 1
                continue

        countline += 1

    return atoms, natom


def connectivity(atoms, natom, bond_max=2.0):
    """
    The algorithm does not consider the connection across periodic boundary. 
    """
    for i in range(natom - 1):
        for j in range(i + 1, natom):
            dist = ((atoms[i].x - atoms[j].x) ** 2 + (atoms[i].y - atoms[j].y) 
                    ** 2 + (atoms[i].z - atoms[j].z) ** 2) ** 0.5
            if dist <= bond_max:
                atoms[i].bonded_atom.append(atoms[j].label)
                atoms[j].bonded_atom.append(atoms[i].label)

    return atoms


def read_data(f1_line, atoms=[], natom=0):
    """
    read Mulliken spin population specific for atoms in supercell. The data set 
    should cover all atomic positions in the supercell. 

    Input is the output file of script Mulliken.py, written in latt_plot.py 
    atomic data format.  
    """
    if natom == 0:
        natom = len(atoms)

    for nline in f1_line[1:]:
        atom_label = int(nline.strip().split()[0])
        atoms[atom_label - 1].data = float(nline.strip().split()[1])

    return atoms


# main function
gui_file = sys.argv[1]
gui = open(gui_file, "r")
gui_line = gui.readlines()
gui.close()

spin_file = sys.argv[2]
spin = open(spin_file, "r")
spin_line = spin.readlines()
spin.close()

atoms, natom = read_geom(gui_line)
atoms = connectivity(atoms, natom)
atoms = read_data(f1_line = spin_line, atoms = atoms, natom = natom)

########################## change according to needs ###########################
neighbor_mx = 11
central_atom = 58
################################################################################
wtspin = open("neighbor_averaged_spin.dat", "w")
wtspin.write('%-8s%12s\n' % ('NEI SEQ', 'MEAN SPIN'))

for i in range(neighbor_mx + 1):
    cluster = Neighbor_cluster(central_atom_lb = central_atom, neighbor_seq = i)
    cluster = cluster.findneighbor(atom_list = atoms)
    # print(cluster.second_lb)
    mean_spin = cluster.average_spin(atom_list = atoms)
    wtspin.write('%-8d%12.8f\n' % (i, mean_spin))

wtspin.close()
