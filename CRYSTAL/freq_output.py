#!/usr/bin/env python
# coding: utf-8

"""
A freq_output class for operating output files of frequency calculations.
Developed only for harmonic phonons.

Required packages: Numpy and os

By Spica. Vir., ICL, May. 21, 22. spica.h.zhou@gmail.com
"""

import os
import numpy as np


class freq_output:
    """
    freq_output class - exact and store key information of frequency output.
    Currently supported attribute:
        self.file: list of string, __init__, Information of output file.
        self.volume: float, get_volume, Lattice volume. Unit: A^3
        self.nmode: int, get_mode_block, Number of vibration modes.
        self.frequency: numpy float array, get_mode, Harmonic vibrational
                        frequency. Unit: THz
        self.eigenvalue: numpy float array, get_mode, Eigenvalues of harmonic
                         modes.Unit: Hartree^2
    """

    def __init__(self, filename):
        """
        Initialize the object.
        Input:
            filename: string, Name of .out file.
        Output:
            self.file
        """
        if os.path.isfile(filename):
            file = open(filename, "r", errors='ignore')
            self.file = file.readlines()
            file.close()
        else:
            print('Error: File not exist.')
            return

    def get_volume(self):
        """
        Obtain lattice information from .out file.
        Input:
            -
        Output:
            self.volume
        """
        label_vol = 'LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - PRIMITIVE CELL'
        countline = 0

        while label_vol not in self.file[countline]:
            countline += 1

        self.volume = float(self.file[countline + 2].strip().split()[-1])

        return

    def get_mode_block(self):
        """
        Get the line indices of frequency information and count the total
        number of available vibration modes.
        Input:
            -
        Output:
            self.nmode
            mode_range, 2*1 numpy int array, line indices of be first and last
                        vibration modes.
        """
        label_bg = 'MODES         EIGV          FREQUENCIES     IRREP  IR   INTENS    RAMAN'
        label_ed = 'NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES (IN BOHR)'
        countline = 0

        while label_bg not in self.file[countline]:
            countline += 1

        range_bg = countline + 2

        while label_ed not in self.file[countline]:
            countline += 1

        range_ed = countline - 2
        nmode = int(self.file[range_ed].strip().split()[1])

        self.nmode = nmode
        mode_range = np.array([range_bg, range_ed], dtype=int)

        return mode_range

    def get_mode(self):
        """
        Get corresponding vibrational frequencies and eigenvalues for all the
        modes.
        Input:
            -
        Output:
            self.frequency
            self.eigenvalue
        """
        if not hasattr(self, 'nmode'):
            mode_range = self.get_mode_block()

        self.frequency = np.array([], dtype=float)
        self.eigenvalue = np.array([], dtype=float)
        for i in self.file[mode_range[0]:mode_range[1] + 1]:
            nm_a = int(i.strip().split()[0].strip('-'))
            nm_b = int(i.strip().split()[1])
            a_freq = float(i.strip().split()[4])
            a_eigv = float(i.strip().split()[2])
            for j in range(nm_a, nm_b + 1):
                self.frequency = np.append(self.frequency, a_freq)
                self.eigenvalue = np.append(self.eigenvalue, a_eigv)

        return
