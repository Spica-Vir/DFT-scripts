#!/usr/bin/env python

################################################################################
# To plot 2D lattice structures stored in .gui file (CRYSTAL fort.34 format).  #
# The color assigned to atoms / bonds can be tuned with a specific value (e.g. #
# population, bondlength, etc.), which is stored in .dat file.                 #
# ----------------------------                                                 #
# by Spica. Vir., ICL, Jun. 21 - 21                                            #
# spica.h.zhou@gmail.com                                                       #
################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Circle
import math
import sys
import copy

class Atom:
    """
    Class atom, plotted as colored / grey circles. Defined by cartesian 
    coordinates. 
        label: int, sequence of atoms listed in .gui file (begin from 1)
        atomic_num: int, atomic number of atoms. Used for determining radius of 
            atoms plotted. Larger atomic_num, larger radius. 
        color: list, assigned by data. Grey as default color. 
        xcoord(ycoord,zcoord): float, cartesian coordinates of atoms. 
        fracx(fracy): float, fractional coordinates (x,y) of atoms. Note, 
            fracx(fracy) correspond to lattice matrix of the input geometry and 
            will be > 1 or < 0 if function cellexpand is used. 
    """

    def __init__(self, label=0, atomic_num=0, 
                 data=None, color=[0.7451, 0.7451, 0.7451], drawn=0, 
                 xcoord=0., ycoord=0., zcoord=0., 
                 fracx=0., fracy=0.):
        from mendeleev import element

        self.label = label
        self.atomic_num = atomic_num
        self.data = data
        self.color = color
        self.drawn = drawn
        self.xcoord = xcoord
        self.ycoord = ycoord
        self.zcoord = zcoord
        self.fracx = fracx
        self.fracy = fracy
        self.atomic_symbol = element(self.atomic_num).symbol

    def drawAtom(self, ax, linecolor='k', linewidth=1, xpos=None, ypos=None):
        """
        xpos(ypos): considering the possible expansion, xpos(ypos) instead of 
            xcoord(ycoord) is used for plotting. Defined by x(y) of Bond object.
            The beginning / ending atoms will be drawn after a bond plotted. 
        """
        if self.atomic_num <= 2:
            radius = 0.2
        elif self.atomic_num > 2 and self.atomic_num <= 10:
            radius = 0.3
        elif self.atomic_num > 10 and self.atomic_num <= 18:
            radius = 0.35
        else:
            radius = 0.40

        if not xpos: 
            xpos = self.xcoord

        if not ypos:
            ypos = self.ycoord

        xy = (xpos, ypos)
        atom = Circle(xy,
                      radius=radius,
                      facecolor=self.color,
                      edgecolor=linecolor,
                      linestyle='-',
                      linewidth=linewidth)
        ax.add_patch(atom)
        atom.set_zorder(1)
        # ax.text(xpos, ypos, self.atomic_symbol, horizontalalignment='center', verticalalignment='center')

class Bond:
    """
    Class bond, plotted as colored / grey lines. Defined by coordinates of the 
    beginning atom and the ending atom. 
        color: list, assigned by data. Grey as default color.
        x: list, [xcoord of beginning atom, ycoord of ending atom]
        y: list, [ycoord of beginning atom, ycoord of ending atom]
        bg(ed)_label: int, label of beginning(ending) atom
        bg(ed)_fx(fy): float, fracx(fracy) of beginning(ending) atom
    """

    def __init__(self, data=None, color=[0.7451, 0.7451, 0.7451], 
                 x=[0., 1.], y=[0., 1.], 
                 bg_label=0, bg_fx=0., bg_fy=0., 
                 ed_label=0, ed_fx=1., ed_fy=1., ):
        self.data = data
        self.color = color
        self.x = x
        self.y = y
        self.bg_label = bg_label
        self.ed_label = ed_label
        self.bg_fx = bg_fx
        self.bg_fy = bg_fy
        self.ed_fx = ed_fx
        self.ed_fy = ed_fy

    def drawBond(self, ax, linewidth=5):
        ax.plot(self.x, self.y, color=self.color, linewidth=linewidth, zorder=0)


class Cell:
    """
    Class Cell, defining the lattice for calculation, i.e., the lattice defined 
    in fort.34 file. Plotted in black dashed line. 
        x(y): 2D-list, the x(y) coordinates of cell corners for plotting. Each 
            element marks the beginning and ending x(y) coordinates of the bond.  
            e.g., x = [[0., x1], [x1, x2], [x2, x3], [x3, x4], [x4, 0.]]
        vect_a(vect_b): list, includes the base 2D lattice vectors. 
            i.e., vect_a = [lattice_matrix[0, 0], lattice_matrix[0, 1]]
                  vect_b = [lattice_matrix[1, 0], lattice_matrix[1, 1]]
    """

    def __init__(self, x=[], y=[], vect_a=[], vect_b=[]):
        self.x = x
        self.y = y
        self.vect_a = vect_a
        self.vect_b = vect_b

    def drawCell(self, ax, linecolor='k', linewidth=0.5, linestyle='dashed'):
        for i in range(len(self.x)):
            ax.plot(self.x[i], self.y[i], color=linecolor,
                    linewidth=linewidth, linestyle=linestyle, zorder = 2)


class Colorbar:
    """
    Class Colorbar, colorbar for plotting atoms / bonds
        ax: figure object, defined by the subplot window
        mx: float, maxima of plotted data
        mi: float, minima of plotted data
        orientation: whether the bar is plotted horizontally or vertically
        colormap: 2D list, each element defines the color(RGB values normalized 
            to 1) at an 'important' point. Colors between those points are 
            obtained by interpolation. 
    """

    def __init__(self, ax, mx=5., mi=0., orientation='vertical', 
                 colormap=[[0, 0, 1], [0, 1, 0],[1, 0, 0]],):
        self.ax = ax
        self.mx = mx
        self.mi = mi
        self.colormap = colormap
        self.orientation = orientation

    def grad_cmap(self, colormap, shrink=10):
        """
        shrink: number of interpolation between 2 neighboring points defined by 
            colormap. 
        """
        npoint = len(colormap)
        new_map = []
        for i in range(1, npoint):
            increment_r = (colormap[i][0] - colormap[i - 1][0]) / shrink
            increment_g = (colormap[i][1] - colormap[i - 1][1]) / shrink
            increment_b = (colormap[i][2] - colormap[i - 1][2]) / shrink
            new_r = colormap[i - 1][0]
            new_g = colormap[i - 1][1]
            new_b = colormap[i - 1][2]
            for j in np.arange(0, 1, 1 / shrink):
                new_r += increment_r
                new_g += increment_g
                new_b += increment_b
                new_map.append([new_r, new_g, new_b])

        new_map.append(colormap[npoint - 1])
        self.colormap = new_map
        return self

    def drawColorbar(self):
        tickpos = np.linspace(self.mi, self.mx, ntick)
        cmap = mpl.colors.ListedColormap(self.colormap)
        norm = mpl.colors.Normalize(vmin=self.mi, vmax=self.mx)
        cb = mpl.colorbar.ColorbarBase(self.ax, cmap=cmap, norm=norm,
                                       orientation=self.orientation)
        cb.set_ticks(tickpos)


def read_geom(gui_line):
    """
    Read the geometry stored in .gui file. Each atom is defined as an object of 
    Atoms class, with geometry data stored. 
        atoms: list, the list of all Atom objects
        natom: int, the number of atoms
        latt_cell: Cell, with 2D lattice vectors and coordinates of corners
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
            cell_x = [[0., float(vect_a[0])], 
                      [float(vect_a[0]), float(vect_a[0]) + float(vect_b[0])],
                      [float(vect_a[0]) + float(vect_b[0]), float(vect_b[0])],
                      [float(vect_b[0]), 0.]]
            cell_y = [[0., float(vect_a[1])],
                      [float(vect_a[1]), float(vect_a[1]) + float(vect_b[1])],
                      [float(vect_a[1]) + float(vect_b[1]), float(vect_b[1])],
                      [float(vect_b[1]), 0.]]
            latt_cell = Cell(x = cell_x, y = cell_y,
                             vect_a = [float(vect_a[0]), float(vect_a[1])],
                             vect_b = [float(vect_b[0]), float(vect_b[1])])
            latvec_matrix = np.matrix([latt_cell.vect_a, latt_cell.vect_b])
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
                an_atom = Atom(label = int(countlabel)%100,
                               atomic_num = int(coord[0]), 
                               xcoord = float(xycoord[0, 0]),
                               ycoord = float(xycoord[0, 1]),
                               zcoord = float(coord[3]),
                               fracx = frac[0, 0], fracy = frac[0, 1])
                atoms.append(an_atom)
            else:
                natom = int(gui_line[countline])
                countline += 1
                continue

        countline += 1

    return atoms, natom, latt_cell

def connectivity(atoms, natom, latt_cell, bond_max=2.0):
    """
    Judge the connectivity between arbitrary 2 atoms, used for building bond
    network. Each bond is defined as a Bond object. 

    Connectivity of bonds across the cell boundary is obtained by building a 2x2
    supercell. See function cellexpand (line 270)

        bond_max: float, maximum allowed distance between 2 atoms regarded to be
            bonded (Angstrom)
        bonds: list, list of all Bond objects
        mean_bond: float, average bond length
    """
    # Bonds cross the cell boundary are plotted on the left & upper boundaries
    # To change that, change parameter expand
    atoms_new, natom_new = cellexpand(latt_cell, atoms = atoms,
                                      expand=[[-1, 1], [0, 2]])
    bonds = []
    mean_bond = 0
    bonded_label = []
    for i in atoms:
        for j in range(natom_new):
            num = atoms_new[j].label
            label_judge = [min([i.label, num]), max([i.label, num])]

            if math.isclose(atoms[num - 1].xcoord, atoms_new[j].xcoord, rel_tol=1e-2) & \
               math.isclose(atoms[num - 1].ycoord, atoms_new[j].ycoord, rel_tol=1e-2):
                rep = 0
            else:
                rep = 1

            if ((rep != 0) & (i.label != num)) | \
            ((rep == 0) & (i.label != num) & (not label_judge in bonded_label)):
                dist = ((i.xcoord - atoms_new[j].xcoord) ** 2 + 
                    (i.ycoord - atoms_new[j].ycoord) ** 2 +
                    (i.zcoord - atoms_new[j].zcoord) ** 2) ** 0.5
                if dist <= bond_max:
                    if rep == 0:
                        bonded_label.append(label_judge)

                    a_bond = Bond(x=[i.xcoord, atoms_new[j].xcoord], 
                                  y=[i.ycoord, atoms_new[j].ycoord], 
                                  bg_label=i.label, 
                                  bg_fx=i.fracx, bg_fy=i.fracy, 
                                  ed_label=atoms_new[j].label,
                                  ed_fx=atoms_new[j].fracx, 
                                  ed_fy=atoms_new[j].fracy)
                    bonds.append(a_bond)
                    mean_bond += dist

    mean_bond /= len(bonds)
    return bonds, mean_bond

def cellexpand(latt_cell, atoms=[], bonds=[], expand=[[0, 1], [0, 1]]):
    """
    Repeating the cell along base cell vectors (2D). Parameter atoms and bonds 
    cannot be simultaneously entered. 

    If atoms is entered, newly introduced atoms will be real objects. Used only 
    for connectivity analysis. If bonds is entered, only data for bonds will be 
    duplicated. Corresponding atoms will be automatically shifted and plotted 
    when plotting. 

        expand: 2D list, defining the beginning and ending points of base 
            lattice vector a and b by fractional numbers. i.e., 
            [[bg_vec_a, ed_vec_a], [bg_vec_b, ed_vec_b]]. The default value is 
            the identical translation. 
        atoms_new(natom_new, bonds_new): expanded atoms(natom, bonds)

    Note: fracx(fracy, bg_fx, bg_fy, ed_fx, ed_fy) is generated based on the 
    original lattice vector, so will be > 1 or < 0. To get the correct cartesian
    coordinates, it should be multiplied by the original lattice vector matrix. 
    """
    # generate lattice points (periodic origins of the original cell)
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
                                       xcoord = xycoord[0, 0], 
                                       ycoord = xycoord[0, 1], zcoord = i.zcoord, 
                                       data = i.data, color = i.color,
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

def read_data(f1_line, f2_line, atoms=[], natom=0, bonds=[]):
    """
    read data specific for atoms/bonds in supercell. The data set needn't to be 
    a complete set for all atoms in supercell. Data are stored in data property 
    of object Atom or Bond correspondingly. By default, atoms/bonds without data
    will be plotted in grey.

    Input data file format: 
    The first line is title line, which will not be read. 

    For atoms data only, input data is stored in variable f1_line. The first
    column should be atom number 1 - N (N is the total number of atoms in 
    supercell). The second column is data for corresponding atoms. 

    For bonds data only, input data is stored in variable f1_line, too. The 1st 
    and 2nd columns are the number of beginning and ending atoms, A1 A2. The 3rd
    column is the data. 

    For atoms + bonds, input data for atoms is stored in variable f1_line, while
    data for bonds is stored in variable f2_line. The formates of data are kept 
    the same as formats above. 
    """
    if natom == 0:
        natom = len(atoms)

    for nline in f1_line[1:]:
        atom_label = int(nline.strip().split()[0])
        atoms[atom_label - 1].data = float(nline.strip().split()[1])

    for nline in f2_line[1:]:
        bg_ed = [int(nline.strip().split()[0]), int(nline.strip().split()[1])]
        bg_label = min(bg_ed)
        ed_label = max(bg_ed)
        for i in bonds:
            bonded_label = [i.bg_label, i.ed_label]
            if [bg_label, ed_label] == [min(bonded_label), max(bonded_label)]:
                i.data = float(nline.strip().split()[2])
            else:
                continue

    return atoms, bonds

def color_assign(data_in, colormap=[[0, 0, 1],
                                    [0, 1, 1],
                                    [0, 1, 0],
                                    [1, 1, 0],
                                    [1, 0, 0]]):
    """
    Defining the color for each atom/bond object acording to its data. 
    For default colormap:
    Red: (1, 0, 0) (max)
    Yellow: (1, 1, 0)
    Green: (0, 1, 0)
    Cyan: (0, 1, 1)
    Blue: (0, 0, 1) (min)

        colormap: 2D list, defined by the list of important points (RGB) and is 
            kept consistent with colormap property of object Colorbar
    """
    nseg = len(colormap) - 1
    data = [i.data for i in data_in if i.data != None]
    if data:
        mx = max(data)
        mi = min(data)
        for i in data_in:
            if i.data != None:
                norm = (i.data - mi) / (mx - mi) * nseg
                pos = float(math.floor(norm))
                if math.isclose(norm, pos, rel_tol=1e-4) and not \
                   math.isclose(pos, 0., rel_tol=1e-4):
                    pos -= 1

                pos = int(pos)
                i.color = [
                    round((norm - pos) * (colormap[pos + 1][0] -
                                    colormap[pos][0]) + colormap[pos][0], 4),
                    round((norm - pos) * (colormap[pos + 1][1] -
                                    colormap[pos][1]) + colormap[pos][1], 4),
                    round((norm - pos) * (colormap[pos + 1][2] -
                                    colormap[pos][2]) + colormap[pos][2], 4)
                ]

            else:
                continue

    return data_in, colormap

def plot(atoms, bonds, mean_bond, latt_cell, out_name='fig_out',
         colormap=[[0, 0, 1], [0, 1, 1], [0, 1, 0], [1, 1, 0], [1, 0, 0]]):
    fig = plt.figure()
    grid = plt.GridSpec(2, 20)
    # plot the bonds and atoms at their two ends. 
    ax = fig.add_subplot(grid[:, 0:19])
    drawn_pos = []
    mx_atomic_num = max([i.atomic_num for i in atoms])
    mi_atomic_num = min([i.atomic_num for i in atoms])
    for i in bonds:
        i.drawBond(ax)
        if not [i.x[0], i.y[0]] in drawn_pos:
            atoms[i.bg_label - 1].drawAtom(ax, xpos=i.x[0], ypos=i.y[0])
            drawn_pos.append([i.x[0], i.y[0]])

        if not [i.x[1], i.y[1]] in drawn_pos:
            atoms[i.ed_label - 1].drawAtom(ax, xpos=i.x[1], ypos=i.y[1])
            drawn_pos.append([i.x[1], i.y[1]])

    # draw colorbar for bonds
    bond_data = [i.data for i in bonds if i.data]
    if bond_data:
        bond_barax = fig.add_subplot(grid[1, 19])
        bond_bar = Colorbar(ax=bond_barax, mx=max(bond_data),
                            mi=min(bond_data))
        bond_bar = bond_bar.grad_cmap(colormap=colormap)
        bond_bar.drawColorbar()

    # draw colorbar for atoms
    atom_data = [i.data for i in atoms if i.data]
    if atom_data:
        atom_barax = fig.add_subplot(grid[0, 19])
        atom_bar = Colorbar(ax=atom_barax, mx=max(atom_data),
                            mi=min(atom_data))
        atom_bar = atom_bar.grad_cmap(colormap=colormap)
        atom_bar.drawColorbar()

    # draw lattice cell
    latt_cell.drawCell(ax)
    # setting up the plot window edges
    rangex = [i[0] for i in drawn_pos]
    rangex.append(max(max(latt_cell.x)))
    rangex.append(min(min(latt_cell.x)))
    rangey = [i[1] for i in drawn_pos]
    rangey.append(max(max(latt_cell.y)))
    rangey.append(min(min(latt_cell.y)))
    xmi = min(rangex) - mean_bond / 4
    xmx = max(rangex) + mean_bond / 4
    ymi = min(rangey) - mean_bond / 4
    ymx = max(rangey) + mean_bond / 4
    ax.set_xlim([xmi, xmx])
    ax.set_ylim([ymi, ymx])
    ax.axis('off')
    ax.set_aspect('equal')

    plt.savefig(fname=out_name + '.pdf', format='pdf')
    plt.show()


# main function
gui_file = sys.argv[1]
gui = open(gui_file, "r")
gui_line = gui.readlines()
gui.close()

file_1 = sys.argv[2]
f1 = open(file_1, "r")
f1_line = f1.readlines()
f1.close()

if len(sys.argv) == 4:
    file_2 = sys.argv[3]
    f2 = open(file_2, "r")
    f2_line = f2.readlines()
    f2.close()
else:
    if len(f1_line[1].strip().split()) == 3:
        f2_line = f1_line
        f1_line = []
    else:
        f2_line = []

atoms, natom, latt_cell = read_geom(gui_line)
bonds, mean_bond = connectivity(atoms, natom, latt_cell)
atoms, bonds = read_data(f1_line, f2_line, atoms, natom, bonds)

atoms, colormap = color_assign(atoms)
bonds, colormap = color_assign(bonds)

# bonds_new = cellexpand(latt_cell, bonds = bonds, expand=[[-1, 1], [-1, 1]])
bonds_new = cellexpand(latt_cell, bonds = bonds, expand=[[0, 1], [0, 1]])

plot(atoms, bonds_new, mean_bond, latt_cell, colormap = colormap)
