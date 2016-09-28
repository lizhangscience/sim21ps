# Module name: SimAGN
# Class: Elp
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
		Params = {'Center':self.Center,'MajAxis':self.MajAxis,\
		'MinAxis':self.MinAxis,'Angle':self.Angle}

	def genCore(self,Rows=512,Cols=512):
		# To avoid outflow, the image mat should be resized based 
		# on the parameters
		AxisMax = max(self.MajAxis,self.MinAxis)
		Row_Add = int(np.ceil(AxisMax))
		Col_Add = int(np.ceil(AxisMax))
		ImageMat = np.zeros((Rows+2*Row_Add,Cols+2*Col_Add))
		
		# Ellipse
		Angs = np.linspace(-np.pi,np.pi,100)
		x = self.MajAxis * np.cos(Angs)
		y = self.MinAxis * np.sin(Angs)
		x_r = self.Center[0] + x*np.cos(self.Angle) - y*np.sin(self.Angle)
		y_r = self.Center[1] + x*np.sin(self.Angle) + y*np.cos(self.Angle)
		# Float to int, to generate the indices
		x_ind = np.round(x_r) - 1
		y_ind = np.round(y_r) - 1
		x_ind = x_ind.astype(np.int)
		y_ind = y_ind.astype(np.int)	
		
		# Fill with the ellipe and resize to given size
		ImageMat[y_ind+Row_Add,x_ind+Col_Add] = 1
		# resize
		ImageCut = ImageMat[Row_Add:Row_Add+Rows,Col_Add:Col_Add+Cols].copy()
				
		return ImageCut

	def genLobes(self,CoreAng = np.pi/2,CoreCen = np.zeros((2,)),Rows=512,Cols=512):
		# Generate two lobes
		# Init
		Angs = np.linspace(-np.pi,np.pi,100)
		# Resize image
		AxisMax = max(self.MajAxis,self.MinAxis)
		Row_Add = int(np.ceil(AxisMax))
		Col_Add = int(np.ceil(AxisMax))
		ImageMat = np.zeros((Rows+2*Row_Add,Cols+2*Col_Add))
		# lobe1
		Rot_Ang = self.Angle + CoreAng
		CenDiff = [self.MajAxis * np.cos(self.Angle),self.MajAxis * np.sin(self.Angle)]
		self.Center[0] = CoreCen[0]+CenDiff[0]*np.cos(CoreAng)-CenDiff[1]*np.sin(CoreAng)
		self.Center[1] = CoreCen[1]+CenDiff[0]*np.sin(CoreAng)+CenDiff[1]*np.cos(CoreAng)
		
		# Gen ellipse
		x_lobe1 = self.MajAxis * np.cos(Angs)
		y_lobe1 = self.MinAxis * np.sin(Angs)
		x_lobe1_r = x_lobe1*np.cos(Rot_Ang)-y_lobe1*np.sin(Rot_Ang)
		y_lobe1_r = x_lobe1*np.sin(Rot_Ang)+y_lobe1*np.cos(Rot_Ang)
		# Transform to indices
		x_lobe1_ind = np.round(x_lobe1_r+self.Center[0]) - 1
		y_lobe1_ind = np.round(y_lobe1_r+self.Center[1]) - 1
		x_lobe1_ind = x_lobe1_ind.astype(np.int)
		y_lobe1_ind = y_lobe1_ind.astype(np.int)

		# love2
		Rot_Ang = self.Angle + CoreAng + np.pi
		CenDiff = [self.MajAxis * np.cos(self.Angle),self.MajAxis * np.sin(self.Angle)]
		self.Center[0] = CoreCen[0]+CenDiff[0]*np.cos(CoreAng+np.pi)-CenDiff[1]*np.sin(CoreAng+np.pi)
		self.Center[1] = CoreCen[1]+CenDiff[0]*np.sin(CoreAng+np.pi)+CenDiff[1]*np.cos(CoreAng+np.pi)
		
		# Gen ellipse
		x_lobe2 = self.MajAxis * np.cos(Angs)
		y_lobe2 = self.MinAxis * np.sin(Angs)
		x_lobe2_r = x_lobe2*np.cos(Rot_Ang)-y_lobe2*np.sin(Rot_Ang)
		y_lobe2_r = x_lobe2*np.sin(Rot_Ang)+y_lobe2*np.cos(Rot_Ang)
		# Transform to indices
		x_lobe2_ind = np.round(x_lobe2_r+self.Center[0]) -1
		y_lobe2_ind = np.round(y_lobe2_r+self.Center[1]) -1
		x_lobe2_ind = x_lobe2_ind.astype(np.int)
		y_lobe2_ind = y_lobe2_ind.astype(np.int)

		# Embed into image
		ImageMat[y_lobe1_ind+Row_Add,x_lobe1_ind+Col_Add] = 1
		ImageMat[y_lobe2_ind+Row_Add,x_lobe2_ind+Col_Add] = 1
		
		# resize image
		ImageCut = ImageMat[Row_Add:Row_Add+Rows,Col_Add:Col_Add+Cols].copy()
				
		return ImageCut 
	
# Calc_Tb
def Calc_Tb(flux_in_Jy,pixel_area,freq):
	c = 2.99792458e8
 	kb = 1.38e-23

	Omegab = pixel_area/(3600*180/np.pi)/(3600*180/np.pi)
	sb = flux_in_Jy * 1e-26 /Omegab
	FluxPixel = Sb/2/freq/freq*c*c/kb

	return FluxPixel

def SimpleSim(Rows=512,Cols=512):
	# Init
	Param_core = Elp()
	Param_lobe = Elp()	
	# Caution: pay attention to the index
	# Core parameters	
	Param_core.Center[0] = np.random.uniform(1,Cols)
	Param_core.Center[1] = np.random.uniform(1,Rows)
	Param_core.MajAxis = np.random.uniform(0,1)
	Param_core.MinAxis = np.random.uniform(0,1)
	Param_core.Angle = np.random.uniform(-np.pi,np.pi)
	# Love parameters
	Param_lobe.MajAxis = np.random.uniform(0,10)
	Param_lobe.MinAxis = np.random.uniform(0,1)
	Param_lobe.Angle = np.random.uniform(-np.pi,np.pi)
	# Embed into the image mat
	ImgLobe = Param_lobe.genLobes(CoreAng=Param_core.Angle, CoreCen=Param_core.Center)
	ImgCore = Param_core.genCore()
	ImageMat = ImgLobe+ImgCore
	# Display
	Idx = np.argwhere(ImageMat>0)
	ImageMat[Idx[:,0],Idx[:,1]] = 100
	ImgTest = Image.fromarray(ImageMat)
	ImgTest.show()

def GenMultiFRs(Rows=512,Cols=512,NumFR=100):
	# Generate multiple simulated FRs
	# Init	
	Param_core = Elp()
	Param_lobe = Elp()	
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
		Param_lobe.MinAxis = np.random.uniform(0,1)
		Param_lobe.Angle = np.random.uniform(-np.pi,np.pi)
		# Embed into the image mat
		ImgLobe = Param_lobe.genLobes(CoreAng=Param_core.Angle, CoreCen=Param_core.Center)
		ImgCore = Param_core.genCore()
		ImageMat = ImageMat + ImgLobe+ImgCore

	# Display
	Idx = np.argwhere(ImageMat>0)
	ImageMat[Idx[:,0],Idx[:,1]] = 100
	ImgTest = Image.fromarray(ImageMat)
	ImgTest.show()








