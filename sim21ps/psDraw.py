#!/usr/bin/env python3
#
# Copyright (c) 2016 Zhixian MA <zxma_sjtu@qq.com>
# MIT license

"""
This module namely psDraw is designed to read csv format ps lists,
calculate the flux and surface brightness of the sources at different
frequency, and then display them on the image mat.

Modules
-------
Basic modules: numpy, pandas, PIL
Custom designed modules: basic_parameters, psCatelogue

Classes
-------
Flux: class
    A class to calculate the ps's surface brightness accordingly
    to its frequency and distance.

Functions
---------
read_csv: read the csv format files, judge the ps type and
transformed to be iterable numpy.ndarray.

calc_flux: calculate the flux and surface brightness of the ps.

draw_elp: processing on the elliptical and circular core or lobes.

draw_ps: draw the ps on the image map.
"""

# Modules
import numpy as np
# from pandas import DataFrame
import pandas as pd
# Cumstom designed modules
# import basic_params
# import psCatelogue

# Init
# Params = basic_params.PixelParams(img_size)


class Flux:
    """
    To calculate the flux and surface brightness of the point sources
    accordingly to its type and frequency

    Parameters
    ----------
    Freq: float
        The frequency
    ClassType: int
        The type of point source, which is default as 1.
        | ClassType | Code |
        |:---------:|:----:|
        |    SF		|  1   |
        |    SB     |  2   |
        |  RQ AGN	|  3   |
        |   FRI     |  4   |
        |   FRII    |  5   |


    Functions
    ---------
    genSpec:
        Generate the spectrum of the source at frequency freq.
    calc_Tb:
        Calculate the average surface brightness, the area of the source
        should be inputed.
    """

    def __init__(self, Freq=150, ClassType=1):

        # Frequency
        self.Freq = Freq
        # ClassType = ClassType
        self.ClassType = ClassType

    def genSpec(self):
        # generate the spectrum
        # Use IF-THEN to replace SWITCH-CASE
        # reference flux at 151MHz, see Willman et al's work
        self.I_151 = 10**(np.random.uniform(-4, -3))
        # Clac flux
        if self.ClassType == 1:
            Spec = (self.Freq / 151e6)**(-0.7) * self.I_151
        elif self.ClassType == 2:
            Spec = (self.Freq / 151e6)**(-0.7) * self.I_151
        elif self.ClassType == 3:
            Spec = (self.Freq / 151e6)**(-0.7) * self.I_151
        elif self.ClassType == 4:
            Spec_lobe = (self.Freq / 151e6)**-0.75 * self.I_151
            a0 = np.log10(self.I_151) - 0.7 * np.log10(151e6) + \
                0.29 * np.log10(151e6) * np.log10(151e6)
            lgs = a0 + 0.7 * np.log10(self.Freq) - 0.29 * \
                np.log10(self.Freq) * np.log10(self.Freq)
            Spec_core = 10**lgs
            Spec = np.array([Spec_core, Spec_lobe])

        elif self.ClassType == 5:
            Spec_lobe = (self.Freq / 151e6)**-0.75 * self.I_151
            Spec_hotspot = (self.Freq / 151e6)**-0.75 * self.I_151
            a0 = np.log10(self.I_151) - 0.7 * np.log10(151e6) + \
                0.29 * np.log10(151e6) * np.log10(151e6)
            lgs = a0 + 0.7 * np.log10(self.Freq) - 0.29 * \
                np.log10(self.Freq) * np.log10(self.Freq)
            Spec_core = 10**lgs
            Spec = np.array([Spec_core, Spec_lobe, Spec_hotspot])

        return Spec

    # calc_Tb
    def calc_Tb(self, Area):
        # light speed
        c = 2.99792458e8
        # ?
        kb = 1.38e-23
        # flux in Jy
        flux_in_Jy = self.genSpec()
        Omegab = Area / (3600 * 180 / np.pi) / (3600 * 180 / np.pi)

        Sb = flux_in_Jy * 1e-26 / Omegab
        FluxPixel = Sb / 2 / self.Freq / self.Freq * c * c / kb

        return FluxPixel


def read_csv(FileName, FoldName='PS_tables'):
    """
    Read csv format point source files,judge its class type according
    to its name.
    For example, 'PS_Num_YYYYMMDD_HHMMSS.csv'
    Split it by '_'

    Parameters
    ----------
    FileName: str
        Name of the file.

    """

    # Split and judge point source type
    ClassList = ['SF', 'SB', 'RQ', 'FRI', 'FRII']
    ClassName = FileName.split('_')[0]
    ClassType = ClassList.index(ClassName) + 1
    # Read csv
    PS_data = pd.read_csv(FoldName + '/' + FileName)

    return ClassType, PS_data


def calc_flux(ClassType, Freq, PS_data):
    """
    Calculate the flux and surface brightness of the point source.

    Parameters
    ----------
    ClassType: int
        Type of point source
    Freq: float
        Frequency
    PS_data: pandas.core.frame.DataFrame
        Data of the point sources
    """
    # init flux
    PS_flux = Flux(Freq=Freq, ClassType=ClassType)
    # PS_flux_list
    NumPS = PS_data.shape[0]
    if ClassType <= 3:
        PS_flux_list = np.zeros((NumPS,))
    else:
        PS_flux_list = np.zeros((NumPS, 2))

    # Iteratively calculate flux
    for i in range(NumPS):
        PS_area = PS_data['Area'][i]
        PS_flux_list[i, :] = PS_flux.calc_Tb(PS_area)[i]

    return PS_flux_list


def draw_rq(ImgMat, PS_data, Freq):
    """
    Designed to draw the radio quiet AGN

    Parameters
    ----------
    ImgMat: np.ndarray
        Two dimensional matrix, to describe the image
    PS_data: pandas.core.frame.DataFrame
        Data of the point sources
    Freq: float
        Frequency
    """

    # Gen flux list
    PS_flux_list = calc_flux(3, Freq, PS_data)
    # Iteratively draw the ps
    NumPS = PS_data.shape[0]
    for i in range(NumPS):
        Core_x = np.ceil(PS_data['Core_x (pix)'][i]) - 1
        Core_y = np.ceil(PS_data['Core_y (pix)'][i]) - 1
        ImgMat[int(Core_y), int(Core_x)] += PS_flux_list[i]

    return ImgMat


def draw_cir(ImgMat, PS_data, ClassType, Freq):
    """
    Designed to draw the circular  star forming  and star bursting PS.

    Prameters
    ---------
    ImgMat: np.ndarray
        Two dimensional matrix, to describe the image
    PS_data: pandas.core.frame.DataFrame
        Data of the point sources
    ClassType: int
        Class type of the point soruces
    Freq: float
        Frequency
    """
    Rows, Cols = ImgMat.shape
    # Gen flux list
    PS_flux_list = calc_flux(ClassType, Freq, PS_data)
    #  Iteratively draw the ps
    NumPS = PS_data.shape[0]
    for i in range(NumPS):
        # grid
        PS_radius = np.ceil(PS_data['radius'][i])  # To be fixed
        Core_x = np.ceil(PS_data['Core_x (pix)'][i]) - 1
        Core_y = np.ceil(PS_data['Core_y (pix)'][i]) - 1
        # Fill with circle
        x = np.arange(Core_x - PS_radius, Core_x + PS_radius + 1, 1)
        y = np.arange(Core_y - PS_radius, Core_y + PS_radius + 1, 1)
        for p in range(len(x)):
            for q in range(len(y)):
                if (x[p] >= 0) and (x[p] <= Cols - 1) and (y[q] >= 0) and (y[q] <= Rows - 1):
                    if np.sqrt((x[p] - Core_x)**2 + (y[q] - Core_y)**2) <= PS_radius:
                        ImgMat[int(y[q]), int(x[p])] += PS_flux_list[i]

    return ImgMat


def draw_lobe(ImgMat, PS_data, ClassType, Freq):
    """
    Designed to draw the elliptical lobes of FRI and FRII

    Prameters
    ---------
    ImgMat: np.ndarray
        Two dimensional matrix, to describe the image
    PS_data: pandas.core.frame.DataFrame
        Data of the point sources
    ClassType: int
        Class type of the point soruces
    Freq: float
        Frequency

    """
    # Init
    Rows, Cols = ImgMat.shape
    NumPS = PS_data.shape[0]
    # Gen flux list
    PS_flux_list = calc_flux(ClassType, Freq, PS_data)

    # Iteratively draw ps
    for i in range(NumPS):
        # Parameters
        Core_x = np.ceil(PS_data['Core_x (pix)'][i]) - 1
        Core_y = np.ceil(PS_data['Core_y (pix)'][i]) - 1
        lobe_maj = np.ceil(PS_data['lobe_maj'][i])
        lobe_min = np.ceil(PS_data['lobe_min'][i])
        lobe_ang = PS_data['lobe_ang'][i]

        # Lobe1
        lobe1_core_x = Core_x + lobe_maj * np.cos(lobe_ang)
        lobe1_core_y = Core_y + lobe_maj * np.sin(lobe_ang)
        # Focuses
        lobe_c = np.sqrt(lobe_maj**2 - lobe_min**2)  # focus distance
        F1_core_x = lobe_c
        F1_core_y = 0
        F2_core_x = -lobe_c
        F2_core_y = 0
        # draw
        a = int(np.round(lobe_maj))
        b = int(np.round(lobe_min))
        x = np.arange(-a, a + 1, 1)
        y = np.arange(-b, b + 1, 1)
        # Ellipse
        for p in range(len(x)):
            for q in range(len(y)):
                DistFocus1 = np.sqrt(
                    (x[p] - F1_core_x)**2 + (y[q] - F1_core_y)**2)
                DistFocus2 = np.sqrt(
                    (x[p] - F2_core_x)**2 + (y[q] - F2_core_y)**2)
                if (DistFocus1 + DistFocus2 <= 2 * lobe_maj):
                    x_r = x[p] * np.cos(lobe_ang) - y[q] * np.sin(lobe_ang)
                    y_r = x[p] * np.sin(lobe_ang) + y[q] * np.cos(lobe_ang)
                    x_r = int(round(x_r + lobe1_core_x))
                    y_r = int(round(y_r + lobe1_core_y))
                    # Judge and Fill
                    if (x_r >= 0) and (x_r <= Cols - 1) and (y_r >= 0) and (y_r <= Rows - 1):
                        ImgMat[int(y_r), int(x_r)] += PS_flux_list[i][1]
        # Lobe2
        lobe2_core_x = Core_x + lobe_maj * np.cos(lobe_ang + np.pi)
        lobe2_core_y = Core_y + lobe_maj * np.sin(lobe_ang + np.pi)
        # Focuses
        lobe_c = np.sqrt(lobe_maj**2 - lobe_min**2)  # focus distance
        F1_core_x = lobe_c
        F1_core_y = 0
        F2_core_x = -lobe_c
        F2_core_y = 0
        # draw
        a = int(np.round(lobe_maj))
        b = int(np.round(lobe_min))
        x = np.arange(-a, a + 1, 1)
        y = np.arange(-b, b + 1, 1)
        # Ellipse
        for p in range(len(x)):
            for q in range(len(y)):
                DistFocus1 = np.sqrt(
                    (x[p] - F1_core_x)**2 + (y[q] - F1_core_y)**2)
                DistFocus2 = np.sqrt(
                    (x[p] - F2_core_x)**2 + (y[q] - F2_core_y)**2)
                if (DistFocus1 + DistFocus2 <= 2 * lobe_maj):
                    x_r = x[p] * np.cos(lobe_ang + np.pi) - \
                        y[q] * np.sin(lobe_ang + np.pi)
                    y_r = x[p] * np.sin(lobe_ang + np.pi) + \
                        y[q] * np.cos(lobe_ang + np.pi)
                    x_r = int(round(x_r + lobe2_core_x))
                    y_r = int(round(y_r + lobe2_core_y))
                    # Judge and Fill
                    if (x_r >= 0) and (x_r <= Cols - 1) and (y_r >= 0) and (y_r <= Rows - 1):
                        ImgMat[int(y_r), int(x_r)] += PS_flux_list[i][1]

        # Core
        ImgMat[int(Core_y), int(Core_x)] += PS_flux_list[i][0]

    return ImgMat


def draw_ps(ImgMat, PS_data, ClassType, Freq):
    """
    Designed to draw the all point sources

    Prameters
    ---------
    ImgMat: np.ndarray
        Two dimensional matrix, to describe the image
    PS_data: pandas.core.frame.DataFrame
        Data of the point sources
    ClassType: int
        Class type of the point soruces
    Freq: float
        Frequency

    """
    if ClassType == 1 or ClassType == 2:
        ImgMat = draw_cir(ImgMat, PS_data, ClassType, Freq)
    elif ClassType == 3:
        ImgMat = draw_rq(ImgMat, PS_data, Freq)
    else:
        ImgMat = draw_lobe(ImgMat, PS_data, ClassType, Freq)

    return ImgMat
