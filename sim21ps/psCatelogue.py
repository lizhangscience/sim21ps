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
import numpy as np
from pandas import DataFrame
# Custom module
import basic_params

# Defination of classes
class Flux:
    def __init__(self,Freq = 150,ClassType=1):
        self.I_151 = 10**(np.random.uniform(-4,-3))
        self.Freq = Freq
        self.ClassType = ClassType

    def genSpec(self):
        # generate the spectrum
        # Use IF-THEN to replace SWITCH-CASE
        if self.ClassType == 1:
            Spec_lobe = (self.Freq/151e6)**-0.75*self.I_151
            a0 = np.log10(self.I_151)-0.7*np.log10(151e6)+0.29*np.log10(151e6)*np.log10(151e6)
            lgs = a0+0.7*np.log10(self.Freq)-0.29*np.log10(self.Freq)*np.log10(self.Freq)
            Spec_core = 10**lgs
            Spec = np.array([Spec_core,Spec_lobe])

        elif self.ClassType == 2:
            Spec_lobe = (self.Freq/151e6)**-0.75*self.I_151
            Spec_hotspot = (self.Freq/151e6)**-0.75*self.I_151
            a0 = np.log10(self.I_151)-0.7*np.log10(151e6)+0.29*np.log10(151e6)*np.log10(151e6)
            lgs = a0+0.7*np.log10(self.Freq)-0.29*np.log10(self.Freq)*np.log10(self.Freq)
            Spec_core = 10**lgs
            Spec = np.array([Spec_core,Spec_lobe,Spec_hotspot])

        return Spec

    # Calc_Tb
    def Calc_Tb(self,Area,Flag=0):
        c = 2.99792458e8
        kb = 1.38e-23
        flux_in_Jy = self.genSpec()[Flag]
        Omegab = Area/(3600*180/np.pi)/(3600*180/np.pi)

        Sb = flux_in_Jy * 1e-26 /Omegab
        FluxPixel = Sb/2/self.Freq/self.Freq*c*c/kb

        return FluxPixel

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
    Lumo: float
        Brightness temperature, [mK] ? 
    
    Functions
    ---------
    gen_sgl_ps
        Generate single ps
    save_as_csv
        
    """
    # Init
    z = 0
    dA = 0
    Columns = []
    nCols = 0
    
    def __init__(self):
        # Redshift
        self.z = np.random.uniform(0,20)
        # angular diameter distance
        Param = basic_params.PixelParams(self.z)
        self.dA = Param.dA
        self.Columns = ['z','dA']
        self.nCols = 2
        
    def gen_sgl_ps(self):
        """
        Generate single point source, and return its data as a list.
        
        """ 
        # Redshift
        self.z = np.random.uniform(0,20)
        # angular diameter distance
        Param = basic_params.PixelParams(self.z)
        self.dA = Param.dA
        PS_list = np.array([self.z,self.dA.value])
        return PS_list
    
    def save_as_csv(self,NumPS = 100):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS,self.nCols))
        for x in range(NumPS):
            PS_Table[x,:] = self.gen_sgl_ps()
        
        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table,columns = self.Columns,index = list(range(NumPS)))
        
        return PS_frame
        
class StarForming(PointSource):
    """
    Generate star forming point sources, inheritate from PointSource class.
    """
    # Init
    z = 0
    dA = 0
    scale = 0 
    Columns = []
    nCols = 0
     
        
    def __init__(self):
       PointSource.__init__ (self)
       self.Columns.append('scale')
       self.nCols += 1
       
    def get_scale(self,Lumo_1400):
       Temp = 0.22 * np.log10(Lumo_1400) - np.log10(1+self.z) - 3.32
       self.scale = 10 ** Temp / 1e3
       
       return self.scale
    
    def gen_sgl_ps(self,Lumo_1400):
        """
        Generate single point source, and return its data as a list.
        
        """ 
        # Redshift
        self.z = np.random.uniform(0,20)
        # angular diameter distance
        Param = basic_params.PixelParams(self.z)
        self.dA = Param.dA
        self.scale = self.get_scale(self,Lumo_1400)
        
        PS_list = np.array([self.z,self.dA.value,self.scale])
        return PS_list
    
    def save_as_csv(self,Lumo_1400,NumPS = 100):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS,self.nCols))
        for x in range(NumPS):
            PS_Table[x,:] = self.gen_sgl_ps(Lumo_1400)
        
        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table,columns = self.Columns,index = list(range(NumPS)))
        
        return PS_frame
        
class StarBursting(PointSource):
    """
    Generate star forming point sources, inheritate from PointSource class.
    """
    # Init
    z = 0
    dA = 0
    scale = 0 
    Columns = []
    nCols = 0
     
        
    def __init__(self):
       PointSource.__init__ (self)
       self.Columns.append('scale')
       self.nCols += 1
       
    def get_scale(self):
        if self.z <= 1.5:
            self.scale = (1 + self.z)**2.5 * 1e-3
        else:
            self.scale = 10 * 1e-3
        
        return self.scale
    
    def gen_sgl_ps(self):
        """
        Generate single point source, and return its data as a list.
        
        """ 
        # Redshift
        self.z = np.random.uniform(0,20)
        # angular diameter distance
        Param = basic_params.PixelParams(self.z)
        self.dA = Param.dA
        self.scale = self.get_scale(self)
        
        PS_list = np.array([self.z,self.dA.value,self.scale])
        return PS_list
    
    def save_as_csv(self,NumPS = 100):
        """
        Generae NumPS of point sources and save them into a csv file.
        """
        # Init
        PS_Table = np.zeros((NumPS,self.nCols))
        for x in range(NumPS):
            PS_Table[x,:] = self.gen_sgl_ps()
        
        # Transform into Dataframe
        PS_frame = DataFrame(PS_Table,columns = self.Columns,index = list(range(NumPS)))
        
        return PS_frame