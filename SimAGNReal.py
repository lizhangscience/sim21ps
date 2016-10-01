# Module name: SimAGN
# Class: Elp
# Class: Flux
# Functions
# 1. Calc_Tb
# 2. SimpleSim
# 3. GenMultiFRs

import numpy as np
import PIL.Image as Image
import pyfits
import matplotlib.pyplot as plt

# define Elliptical lobe and core class
class Elp:
    def __init__(self):
        self.Center = np.zeros((2,))
        self.MajAxis = np.zeros((1,))
        self.MinAxis = np.zeros((1,))
        self.Angle = np.zeros((1,))

    def save_as_dict(self):
        Params = {'Center':self.Center,'MajAxis':self.MajAxis,
                  'MinAxis':self.MinAxis,'Angle':self.Angle}
        return Params

    def genCore(self,ImageMat,Param_Flux):
        # preparing
        Area = np.pi * self.MajAxis * self.MinAxis
        Rows,Cols = ImageMat.shape
        if self.MajAxis > self.MinAxis:
            AxisMax = self.MajAxis
            C_core = np.sqrt(self.MajAxis**2 - self.MinAxis**2)
            # Focus1
            F1_core_x = C_core
            F1_core_y = 0
            # Focus2
            F2_core_x = -C_core
            F2_core_y = 0
        else:
            AxisMax = self.MinAxis
            C_core = np.sqrt(self.MinAxis**2 - self.MajAxis**2)
            # Focus1
            F1_core_x = 0
            F1_core_y = C_core
            # Focus2
            F2_core_x = 0
            F2_core_y = -C_core

        # Fill with flux
        # Intergerize
        a = int(np.round(self.MajAxis))
        b = int(np.round(self.MinAxis))
        x = np.arange(-a,a+1,1)
        y = np.arange(-b,b+1,1)
        # Ellipse
        for i in range(len(x)):
            for j in range(len(y)):
                DistFocus1 = np.sqrt((x[i]-F1_core_x)**2+(y[j]-F1_core_y)**2)
                DistFocus2 = np.sqrt((x[i]-F2_core_x)**2+(y[j]-F2_core_y)**2)
                if     (DistFocus1+DistFocus2<=2*AxisMax):
                    x_r = x[i]*np.cos(self.Angle) - y[j]*np.sin(self.Angle)
                    y_r = x[i]*np.sin(self.Angle) + y[j]*np.cos(self.Angle)
                    x_r = int(round(x_r+self.Center[0]))
                    y_r = int(round(y_r+self.Center[1]))
                    # Judge and Fill
                    if (x_r>=1) and (x_r<=Cols) and (y_r>=1) and (y_r<=Rows):
                        ImageMat[y_r-1][x_r-1] = Param_Flux.Calc_Tb(Area,Flag=0)

        return ImageMat

    def genLobes(self,ImageMat,Param_Flux,CoreAng=np.pi/2,CoreCen=np.zeros((2,))):
        # preparing
        Area = np.pi * self.MajAxis * self.MinAxis
        Rows,Cols = ImageMat.shape
        # Lobe1
        RotAng = self.Angle + CoreAng
        CenDiff = [self.MajAxis * np.cos(self.Angle),self.MajAxis * np.sin(self.Angle)]
        self.Center[0] = CoreCen[0]+CenDiff[0]*np.cos(CoreAng)-CenDiff[1]*np.sin(CoreAng)
        self.Center[1] = CoreCen[1]+CenDiff[0]*np.sin(CoreAng)+CenDiff[1]*np.cos(CoreAng)

        if self.MajAxis > self.MinAxis:
            AxisMax = self.MajAxis
            C_core = np.sqrt(self.MajAxis**2 - self.MinAxis**2)
            # Focus1
            F1_core_x = C_core
            F1_core_y = 0
            # Focus2
            F2_core_x = -C_core
            F2_core_y = 0
        else:
            AxisMax = self.MinAxis
            C_core = np.sqrt(self.MinAxis**2 - self.MajAxis**2)
            # Focus1
            F1_core_x = 0
            F1_core_y = C_core
            # Focus2
            F2_core_x = 0
            F2_core_y = -C_core

        a = int(np.round(self.MajAxis))
        b = int(np.round(self.MinAxis))
        x = np.arange(-a,a+1,1)
        y = np.arange(-b,b+1,1)
        # Ellipse
        for i in range(len(x)):
            for j in range(len(y)):
                DistFocus1 = np.sqrt((x[i]-F1_core_x)**2+(y[j]-F1_core_y)**2)
                DistFocus2 = np.sqrt((x[i]-F2_core_x)**2+(y[j]-F2_core_y)**2)
                if     (DistFocus1+DistFocus2<=2*AxisMax):
                    x_r = x[i]*np.cos(RotAng) - y[j]*np.sin(RotAng)
                    y_r = x[i]*np.sin(RotAng) + y[j]*np.cos(RotAng)
                    x_r = int(round(x_r+self.Center[0]))
                    y_r = int(round(y_r+self.Center[1]))
                    # Judge and Fill
                    if (x_r>=1) and (x_r<=Cols) and (y_r>=1) and (y_r<=Rows):
                        ImageMat[y_r-1][x_r-1] = Param_Flux.Calc_Tb(Area,Flag=1)

        # Lobe2
        Rot_Ang = self.Angle + CoreAng + np.pi
        CenDiff = [self.MajAxis * np.cos(self.Angle),self.MajAxis * np.sin(self.Angle)]
        self.Center[0] = CoreCen[0]+CenDiff[0]*np.cos(CoreAng + np.pi)-CenDiff[1]*np.sin(CoreAng + np.pi)
        self.Center[1] = CoreCen[1]+CenDiff[0]*np.sin(CoreAng + np.pi)+CenDiff[1]*np.cos(CoreAng + np.pi)

        if self.MajAxis > self.MinAxis:
            AxisMax = self.MajAxis
            C_core = np.sqrt(self.MajAxis**2 - self.MinAxis**2)
            # Focus1
            F1_core_x = C_core
            F1_core_y = 0
            # Focus2
            F2_core_x = -C_core
            F2_core_y = 0
        else:
            AxisMax = self.MinAxis
            C_core = np.sqrt(self.MinAxis**2 - self.MajAxis**2)
            # Focus1
            F1_core_x = 0
            F1_core_y = C_core
            # Focus2
            F2_core_x = 0
            F2_core_y = -C_core

        a = int(np.round(self.MajAxis))
        b = int(np.round(self.MinAxis))
        x = np.arange(-a,a+1,1)
        y = np.arange(-b,b+1,1)
        # Ellipse
        for i in range(len(x)):
            for j in range(len(y)):
                DistFocus1 = np.sqrt((x[i]-F1_core_x)**2+(y[j]-F1_core_y)**2)
                DistFocus2 = np.sqrt((x[i]-F2_core_x)**2+(y[j]-F2_core_y)**2)
                if     (DistFocus1+DistFocus2<=2*AxisMax):
                    x_r = x[i]*np.cos(RotAng) - y[j]*np.sin(RotAng)
                    y_r = x[i]*np.sin(RotAng) + y[j]*np.cos(RotAng)
                    x_r = int(round(x_r+self.Center[0]))
                    y_r = int(round(y_r+self.Center[1]))
                    # Judge and Fill
                    if (x_r>=1) and (x_r<=Cols) and (y_r>=1) and (y_r<=Rows):
                        ImageMat[y_r-1][x_r-1] = Param_Flux.Calc_Tb(Area,Flag=1)

        return ImageMat


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

def SimpleSim(Rows=512,Cols=512):
    # Init
    Param_core = Elp()
    Param_lobe = Elp()
    Param_Flux = Flux()
    ImageMat = np.zeros((Rows,Cols))
    # Caution: pay attention to the index
    # Core parameters
    Param_core.Center[0] = np.random.uniform(1,Cols)
    Param_core.Center[1] = np.random.uniform(1,Rows)
    Param_core.MajAxis = np.random.uniform(0,1)
    Param_core.MinAxis = np.random.uniform(0,1)
    Param_core.Angle = np.random.uniform(-np.pi,np.pi)
    # Love parameters
    Param_lobe.MajAxis = np.random.uniform(0,10)
    Param_lobe.MinAxis = np.random.uniform(0,4)
    Param_lobe.Angle = np.random.uniform(-np.pi,np.pi)
    # Embed into the image mat
    ImgLobe = Param_lobe.genLobes(ImageMat,Param_Flux,CoreAng=Param_core.Angle, CoreCen=Param_core.Center)
    ImgCore = Param_core.genCore(ImageMat,Param_Flux)
    ImageMat = ImgLobe+ImgCore
    # Display
    #Idx = np.argwhere(ImageMat>0)
    #ImageMat[Idx[:,0],Idx[:,1]] = 100
    ImgTest = Image.fromarray(ImageMat)
    ImgTest.show()

def GenMultiFRs(Rows=512,Cols=512,Freq=150,NumFR=100):
    # Generate multiple simulated FRs
    # Init
    Param_core = Elp()
    Param_lobe = Elp()
    Param_Flux = Flux()
    Param_Flux.Freq = Freq
    ImageMat = np.zeros((Rows,Cols))

    for x in range(NumFR):
        print 'FR %d' % x
        # Core parameters
        Param_core.Center[0] = np.random.uniform(1,Cols)
        Param_core.Center[1] = np.random.uniform(1,Rows)
        Param_core.MajAxis = np.random.uniform(0,1)
        Param_core.MinAxis = np.random.uniform(0,1)
        Param_core.Angle = np.random.uniform(-np.pi,np.pi)
        # Lobe parameters
        Param_lobe.MajAxis = np.random.uniform(0,5)
        Param_lobe.MinAxis = np.random.uniform(0,2)
        Param_lobe.Angle = np.random.uniform(-np.pi,np.pi)
        # Embed into the image mat
        ImgLobe = Param_lobe.genLobes(ImageMat,Param_Flux,CoreAng=Param_core.Angle, CoreCen=Param_core.Center)
        ImgCore = Param_core.genCore(ImageMat,Param_Flux)
        ImageMat = ImgLobe+ImgCore

    # Display
    Idx = np.argwhere(ImageMat>0)
    ImageMat[Idx[:,0],Idx[:,1]] = 100
    ImgTest = Image.fromarray(ImageMat)
    ImgTest = ImgTest.convert('RGB')
    FileName = 'Img_'+str(Freq)+'.jpg'
    ImgTest.save(FileName)
    ImgTest.show()








