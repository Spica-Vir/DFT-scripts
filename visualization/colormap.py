#!/usr/bin/env python
# coding: utf-8

# A Colorscale class to assist plotting a pseudo color figure.
# By Spica. Vir., ICL, Apr. 17, 22. spica.h.zhou@gmail.com

import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.colors import Normalize
from matplotlib.colorbar import Colorbar


class Colorscale:
    """
    Define a Colorscale class to create a pseudo color figure. Methods:
        1. Generate a user-defined colormap
        2. Linearly scale the input data to generate a corresponding color
        3. Assign grey color to empty input data/equal bottom & top values
        4. Export the defined colormap as a Colormap object in matplotlib
        5. Generate a colorbar object according to the user-defined colormap
    By Spica. Vir., ICL, Apr. 17, 22. spica.h.zhou@gmail.com
    """

    def __init__(self, divid=100, colormap=[[0, 0, 1],
                                            [0, 1, 1],
                                            [0, 1, 0],
                                            [1, 1, 0],
                                            [1, 0, 0]], botm=0, top=99):
        """
        divid: int, the accuracy of interpolation, i.e., the maximum number of
               colors plotted in a figure. Divid is recommended to be the
               integer multiply of the number of rows of the colormap matrix.
        colormap: list, define the key points of a user defined colormap. The
                  default map (bottom to top): blue-cyan-green-yellow-red
        botm: float, the value at the bottom of colorbar
        top: float, the value on the top of colorbar
        """
        colormap = np.matrix(colormap, dtype=float)
        colormap = colormap.transpose()

        if divid < np.size(colormap, 1) - 1:
            print('Warning: The number of interpolation is smaller than the \
                   size of colormap. Colormap size will be used.')
            divid = np.size(colormap, 1) - 1

        div_per_range = int(divid // (np.size(colormap, 1) - 1))

        colormx = []
        for i in range(3):
            colorarray = []
            for j in range(np.size(colormap, 1) - 1):
                rangearray = np.linspace(colormap[i, j],
                                         colormap[i, j + 1],
                                         div_per_range + 1)
                rangearray = rangearray[: -1].tolist()
                colorarray = colorarray + rangearray

            colorarray.append(colormap[i, -1])
            colormx.append(colorarray)

        colormx = np.around(np.matrix(colormx, dtype=float), 3)

        self.colormap = colormx.transpose().tolist()
        self.top = top
        self.botm = botm
        self.divid = divid

    def scale_color(self, inpdata):
        """
        Assign the color to input data.
        Data is linearly scaled. Grey is assigned to data exceeding the range.
        """
        if self.botm == self.top:
            data_color = [0.7451, 0.7451, 0.7451]
            return data_color

        color_pos = round((inpdata - self.botm) / (self.top - self.botm)
                          * self.divid)

        if color_pos not in range(len(self.colormap)):
            print('Warning: Data: ', inpdata, ' exceeding pre-defined range. \
                   Grey will be assigned')
            data_color = [0.7451, 0.7451, 0.7451]
        else:
            data_color = self.colormap[color_pos]

        return data_color

    def exportcmap(self, mapname='user_map'):
        """
        Export colormap object with a user-defined name.
        """
        cmap = ListedColormap(self.colormap, name=mapname)

        return cmap

    def draw_colorbar(self, plotax, title=None,
                      ntick=5, ticklabel=[], orientation='Vertical'):
        """
        Draw colorbar on a figure.
        plotax: axis, the axis for plotting.
        title: str, title / label of colorbar
        ntick: list, number of ticks. Top & bottom of a colorbar are marked.
        ticklabel: list, label of ticks.
        orientation: orientation of colorbar, 'Vertical' / 'Horizontal'

        """
        if self.botm == self.top:
            print('Warning: No variations in input data.')
            return

        tickpos = np.linspace(self.botm, self.top, ntick)

        if ticklabel and len(ticklabel) != ntick:
            print('Warning: The dimension of list \'ticklabel\' should equal \
                   ntick. Values are used instead.')
            ticklabel = np.around(tickpos, 2).tolist()
        elif not ticklabel:
            ticklabel = np.around(tickpos, 2).tolist()

        cmap = ListedColormap(self.colormap)
        norm = Normalize(vmin=self.botm, vmax=self.top)
        colorbar = Colorbar(plotax, cmap=cmap, norm=norm,
                            orientation=orientation)
        colorbar.set_ticks(tickpos, labels=ticklabel)

        if title:
            colorbar.set_label(title)

        return colorbar
