#!/usr/bin/env python3
#
# Copyright (c) 2016 Zhixian MA <zxma_sjtu@qq.com>
# MIT license

"""
Simulation of point sources for 21cm signal detection

Point sources types
-----
1. Star forming (SF) galaxies
2. Star bursting (SB) galaxies
3. Radio quiet AGN (RQ_AGN)
4. Faranoff-Riley I (FRI)
5. Faranoff-Riley II (FRII)

References
------------
[1] Wilman R J, Miller L, Jarvis M J, et al.
    A semi-empirical simulation of the extragalactic radio continuum sky for
    next generation radio telescopes[J]. Monthly Notices of the Royal
    Astronomical Society, 2008, 388(3):1335–1348.
[2] Jelic et al.

"""
# Import packages or modules
import os
import numpy as np
import time
from pandas import DataFrame
# Custom module
import basic_params

# Defination of classes
class PointSource:
    """
    The basic class of point sources

    Parameters
    ----------
    z: float
        Redshift, z ~ U(0,20)
    dA: float
        Angular diameter distance, which is calculated according to the cosmology
        constants. In this work, it is calculated by module basic_params
    core_x,core_y: float
        Core coordinates of the point source, which obey to uniform distribution.
    Lumo: float
        Brightness temperature, [mK] ?
        area: float
                Area of the point sources, pix^2

    Functions
    ---------
    gen_sgl_ps
        Generate single ps
    save_as_csv

    """
    # Init
    img_size = 0
    z = 0
    dA = 0
    core_x = 0
    core_y = 0
    area = 1
    Columns = []
    nCols = 0

    def __init__(self, img_size=512, ang_total=1):
        # Image size
        self.img_size = img_size
        # Redshift
        self.z = np.random.uniform(0, 20)
        # angular diameter distance
        self.ang_total = ang_total
        self.Param = basic_params.PixelParams( self.z, self.img_size, self.ang_total)
        self.dA = self.Param.dA
        # Area
        self.area = 1
        # PS_list information
        self.Columns = ['z', 'dA (Mpc)', 'Core_x (pix)',
                        'Core_y (pix)', 'Area']
        self.nCols = 5

    def gen_sgl_ps(self):
        """
        Generate single point source, and return its data as a list.

        """
        # Redshift
        self.z = np.random.uniform(0, 20)
        # angular diameter distance
        self.Param = basic_params.PixelParams(
            self.z, self.img_size, self.ang_total)
        self.dA = self.Param.dA
        # Position
        self.core_x = np.random.uniform(self.img_size - 1) + 1
        self.core_y = np.random.uniform(self.img_size - 1) + 1

        PS_list = np.array(
            [self.z, self.dA.value, self.core_x, self.core_y, self.area])
        return PS_list

    def save_as_csv(self, NumPS=100, folder_name='PS_tables/'):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS, self.nCols))
        for x in range(NumPS):
            PS_Table[x, :] = self.gen_sgl_ps()

        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table, columns=self.Columns,
                             index=list(range(NumPS)))

        # Save to csv
        if os.path.exists(folder_name) == False:
            os.mkdir(folder_name)

        file_name = 'PS_' + str(NumPS) + '_' + \
            time.strftime('%Y%m%d_%H%M%S') + '.csv'
        PS_frame.to_csv(folder_name + '/' + file_name)
        return PS_frame


class StarForming(PointSource):
    """
    Generate star forming point sources, inheritate from PointSource class.
    """
    # Init
    radius = 0
    Lumo_1400 = 0

    def __init__(self, Lumo_1400=1500, img_size=512, ang_total=1):
        PointSource.__init__(self, img_size, ang_total)
        self.Lumo_1400 = Lumo_1400
        self.Columns.append('radius')
        self.nCols += 1

    def get_radius(self):
        Temp = 0.22 * np.log10(self.Lumo_1400) - np.log10(1 + self.z) - 3.32
        self.radius = 10 ** Temp / 2

        return self.scale

    def gen_sgl_ps(self):
        """
        Generate single point source, and return its data as a list.

        """
        # Redshift
        self.z = np.random.uniform(0, 20)
        # angular diameter distance
        self.Param = basic_params.PixelParams(
        self.z, self.img_size, self.ang_total)
        self.dA = self.Param.dA
        self.radius = self.Param.get_angle(self.get_radius())[0]
        # Area
        self.area = np.pi * self.radius**2
        # Position
        self.core_x = np.random.uniform(self.img_size - 1) + 1
        self.core_y = np.random.uniform(self.img_size - 1) + 1

        PS_list = np.array(
            [self.z, self.dA.value, self.core_x, self.core_y, self.area, self.radius])
        return PS_list

    def save_as_csv(self, NumPS=100, folder_name='PS_tables'):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS, self.nCols))
        for x in range(NumPS):
            PS_Table[x, :] = self.gen_sgl_ps()

        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table, columns=self.Columns,
                             index=list(range(NumPS)))

        # Save to csv
        if os.path.exists(folder_name) == False:
            os.mkdir(folder_name)

        file_name = 'SF_' + str(NumPS) + '_' + \
            time.strftime('%Y%m%d_%H%M%S') + '.csv'
        PS_frame.to_csv(folder_name + '/' + file_name)
        return PS_frame


class StarBursting(PointSource):
    """
    Generate star forming point sources, inheritate from PointSource class.
    """
    # Init
    radius = 0

    def __init__(self, img_size=512, ang_total=1):
        PointSource.__init__(self, img_size, ang_total)
        self.Columns.append('radius')
        self.nCols += 1

    def get_radius(self):
        if self.z <= 1.5:
            self.radius = (1 + self.z)**2.5 * 1e-3
        else:
            self.radius = 10 * 1e-3

        return self.radius

    def gen_sgl_ps(self):
        """
        Generate single point source, and return its data as a list.

        """
        # Redshift
        self.z = np.random.uniform(0, 20)
        # angular diameter distance
        self.Param = basic_params.PixelParams(
            self.z, self.img_size, self.ang_total)
        self.dA = self.Param.dA
        self.radius = self.Param.get_angle(self.get_radius())[0]
        # Area
        self.area = np.pi * self.radius**2
        # Position
        self.core_x = np.random.uniform(self.img_size - 1) + 1
        self.core_y = np.random.uniform(self.img_size - 1) + 1

        PS_list = np.array(
            [self.z, self.dA.value, self.core_x, self.core_y, self.area, self.radius])
        return PS_list

    def save_as_csv(self, NumPS=100, folder_name='PS_tables'):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS, self.nCols))
        for x in range(NumPS):
            PS_Table[x, :] = self.gen_sgl_ps()

        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table, columns=self.Columns,
                             index=list(range(NumPS)))

        # Save to csv
        if os.path.exists(folder_name) == False:
            os.mkdir(folder_name)

        file_name = 'SB_' + str(NumPS) + '_' + \
            time.strftime('%Y%m%d_%H%M%S') + '.csv'
        PS_frame.to_csv(folder_name + '/' + file_name)

        return PS_frame


class FRI(PointSource):
    """
    Generate Faranoff-Riley I (FRI) AGN

    Parameters
    ----------
    lobe_maj: float
        The major half axis of the lobe
    lobe_min: float
        The minor half axis of the lobe
    lobe_ang: float
        The rotation angle of the lobe from LOS

    """
    # New parameters
    lobe_maj = 0
    lobe_min = 0
    lobe_ang = 0

    def __init__(self):
        PointSource.__init__(self, img_size=512, ang_total=1)
        self.Columns.extend(
            ['lobe_maj', 'lobe_min', 'lobe_ang'])
        self.nCols += 3

    def gen_lobe(self):
        """
        According to Wang's work, the linear scale at redshift z obeys to U(0,D0(1+z)^(-1.4))
        """
        D0 = 1
        self.lobe_maj = 0.5 * np.random.uniform(0, D0 * (1 + self.z)**(-1.4))
        self.lobe_min = self.lobe_maj * np.random.uniform(0.2, 1)
        self.lobe_ang = np.random.uniform(0, np.pi)

        # Transform to pixel
        self.lobe_maj = self.Param.get_angle(self.lobe_maj)[0]
        self.lobe_min = self.Param.get_angle(self.lobe_min)[0]
        lobe = [self.lobe_maj, self.lobe_min, self.lobe_ang]

        return lobe

    def gen_sgl_ps(self):
        """
        Generate single point source, and return its data as a list.

        """
        # Redshift
        self.z = np.random.uniform(0, 20)
        # angular diameter distance
        self.Param = basic_params.PixelParams(self.z, self.img_size, self.ang_total)
        self.dA = self.Param.dA

        # Position
        self.core_x = np.random.uniform(self.img_size - 1) + 1
        self.core_y = np.random.uniform(self.img_size - 1) + 1

        # Area
        self.area = np.pi * self.lobe_maj * self.lobe_min

        # lobe
        lobe = self.gen_lobe()

        PS_list = [self.z, self.dA.value, self.core_x, self.core_y, self.area]
        PS_list.extend(lobe)

        PS_list = np.array(PS_list)
        return PS_list

    def save_as_csv(self, NumPS=100, folder_name='PS_tables'):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS, self.nCols))
        for x in range(NumPS):
            PS_Table[x, :] = self.gen_sgl_ps()

        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table, columns=self.Columns,
                             index=list(range(NumPS)))

        # Save to csv
        if os.path.exists(folder_name) == False:
            os.mkdir(folder_name)

        file_name = 'FRI_' + str(NumPS) + '_' + \
            time.strftime('%Y%m%d_%H%M%S') + '.csv'
        PS_frame.to_csv(folder_name + '/' + file_name)

        return PS_frame


class FRII(FRI):
    """
    Generate Faranoff-Riley I (FRI) AGN, a class inherit from FRI
    """

    def __init__(self,img_size = 512,ang_total = 1):
        FRI.__init__(self,img_size,ang_total)

    def gen_lobe(self):
        """
        According to Wang's work, the linear scale at redshift z obeys to U(0,D0(1+z)^(-1.4))
        """
        D0 = 1
        self.lobe_maj = 0.5 * np.random.uniform(0, D0 * (1 + self.z)**(-1.4))
        self.lobe_min = self.lobe_maj * np.random.uniform(0.2, 1)
        self.lobe_ang = np.random.uniform(0, np.pi / 3)  # Different from FRI

        # Transform to pixel
        self.lobe_maj = self.Param.get_angle(self.lobe_maj)[0]
        self.lobe_min = self.Param.get_angle(self.lobe_min)[0]
        lobe = [self.lobe_maj, self.lobe_min, self.lobe_ang]

        return lobe
