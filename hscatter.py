#!/bin/python
import sgems
import matplotlib.pyplot as plt
from scipy import stats
from matplotlib import colors,ticker,cm
from scipy.interpolate import griddata
from scipy import interpolate
import numpy as np
import math
import random

#...........................................................
def retire_nan (x,y,z,v):
	xn =[]
	yn =[]
	zn =[]
	vn =[]
	for i in range(0,len(v)):
		if (v[i] > -99999999999999999):
			
			xn.append(x[i])
			yn.append(y[i])
			zn.append(z[i])
			vn.append(v[i])
	return xn, yn, zn, vn
#...........................................................			


"""
VARIOMAP
"""

class   hscatter:
    def _init_(self):
        pass 
    def initialize(self,params):
        self.params = params
        return True
    def execute(self):

	'''
		HSCATTERPLOT 
		_____________________________________________________________________________

		This is a template for SGems program variogram and covariance maps.
		The user have to fill the parameters with a ui.file with same name of this python 
		file before start program. The QT userinterface connect SGems data with this routine
		using the class params.
		
		The output of this program is a graphic demonstranting interpolated variogram values 
		along one plane 

	AUTHOR: DAVID ALVARENGA DRUMOND 		2016
	
	'''	

	'''
	INITIAL PARAMETERS 
	______________________________________________________________________________

	All initial parameters are transcribe by ui.file interface. Variables are choose 
	among common variables. The follow variables are:

	COMMON VARIABLES: 

	head_property= head property 
	head_grid = head grid 
	tail_property = tail property 
	tail_grid= tail grid 
	C_variogram = check box for variogram maps 
	C_covariance= check box for covariance maps 
	nlags= number of lags in experimental variogram 
	lagdistance = lag length for this experimental variograms 
	lineartolerance= linear tolerance for this map
	htolerance = horizontal tolerance for this map 
	vtolerance = vertical tolerance for this lag
	hband = horizontal band 
	vband = vertical band 
	dip = dip of the rake of variogram map 
	azimute = azimuth of the rake of variogram map
	gdiscrete = number of cells to discretize map 
	ncontour = number of contourn lines to interpolate in map 

	
	'''


	# Stabilish initial parameters 

	head_property= self.params['head_prop']['property']
	head_grid = self.params['head_prop']['grid']
	tail_property = self.params['tail_prop']['property']
	tail_grid= self.params['tail_prop']['grid']
	nlags= int(self.params['nlags']['value'])
	lagdistance = float(self.params['lagdistance']['value'])
	lineartolerance= float(self.params['lineartolerance']['value'])
	htolerance = float(self.params['htolangular']['value'])
	vtolerance = float(self.params['vtolangular']['value'])
	hband = float(self.params['hband']['value'])
	vband = float(self.params['vband']['value'])
	Dip = float(self.params['Dip']['value'])
	Azimute= float(self.params['Azimute']['value'])

	print (self.params)
	
	# Get values from SGems and exclude nan values  
	
	xh = []
	yh = []
	zh = [] 
	vh = []
	
	xh1 = sgems.get_property(head_grid, "_X_")
	yh1 = sgems.get_property(head_grid, "_Y_")
	zh1 = sgems.get_property(head_grid, "_Z_")
	vh1 = sgems.get_property(head_grid, head_property)	

	xh, yh, zh, vh = retire_nan(xh1,yh1,zh1,vh1)


	xt = []
	yt = []
	zt = []	
	vt = []

	xt1 = sgems.get_property(tail_grid, "_X_")
	yt1 = sgems.get_property(tail_grid, "_Y_")
	zt1 = sgems.get_property(tail_grid, "_Z_")
	vt1 = sgems.get_property(tail_grid, tail_property)

	xt, yt, zt, vt = retire_nan(xt1,yt1,zt1,vt1)	

	# Verify number of dimensions of problem
	
	# Verify dimensions of head 

	cdimensaoxh = False
	cdimensaoyh = False
	cdimensaozh = False

	for x in range(0,len(xh)):
		if (xh[x] != 0):
			cdimensaoxh = True


	for y in range(0,len(yh)):
		if (yh[y] != 0):
			cdimensaoyh = True


	for z in range(0,len(zh)):
		if (zh[z] != 0):
			cdimensaozh = True

	# Verify dimensions of tail 

	cdimensaoxt = False
	cdimensaoyt = False
	cdimensaozt = False

	for x in range(0,len(xt)):
		if (xt[x] != 0):
			cdimensaoxt = True

	for y in range(0,len(yt)):
		if (yt[y] != 0):
			cdimensaoyt = True


	for z in range(0,len(zt)):
		if (zt[z] != 0):
			cdimensaozt = True

	if ((cdimensaoxt == cdimensaoxh) and (cdimensaoyt == cdimensaoyh) and (cdimensaozt == cdimensaozh)):
		


		# Define maximum distance permissible 
		
		dmaximo = (nlags + 0.5)*lagdistance
		
		# Define tolerances if they are zero 
 
		if (htolerance == 0):
			htolerance= 45
		if (vtolerance == 0):
			vtolerance = 45
		if (hband == 0):
			hband = lagdistance/2
		if (vband == 0):
			vband = lagdistance/2
		if (lagdistance == 0):		
			lagdistance = 100
		if (lineartolerance == 0):
			lineartolerance = lagdistance/2
		
		
		
	
		# Convert tolerances in radians 

		htolerance = math.radians(htolerance)
		vtolerance = math.radians(vtolerance)	

		# Determine dip and azimuth projections  
		htolerance = math.cos(htolerance)
		vtolerance = math.cos(vtolerance)
	
		# Define all euclidian distances permissible 
		
		cabeca = []
		rabo = []
		distanciah =[]	
		distanciaxy =[]	
		distanciax =[]
		distanciay = []
 		distanciaz = []
		for t in range(0, len(yt)):
			for h in range(t, len(yh)):
				cabeca.append(vh[h])
				rabo.append(vt[t])
				dx = xh[h]-xt[t] 				
				dy= yh[h]-yt[t]
				dz = zh[h]-zt[t]
				if (distanciah > 0):
					distanciay.append(dy)
					distanciax.append(dx)
					distanciaz.append(dz)
					distanciah.append(math.sqrt(math.pow(dy,2) + math.pow(dx,2) + math.pow(dz,2)))						
					distanciaxy.append(math.sqrt(math.pow(dy,2) + math.pow(dx,2)))

		# Calculate all cos and sin values 		
	
		cos_Azimute = math.cos(math.radians(90-Azimute))
		sin_Azimute = math.sin(math.radians(90-Azimute))
		sin_Dip = math.sin(math.radians(90-Dip))
		cos_Dip = math.cos(math.radians(90-Dip))
	
		v_valores_admissiveis_h =[]
		v_valores_admissiveis_t =[]
		

		# Calculate admissible points  
		
		

		for l in range(0, nlags):			
			valores_admissiveis_h = []
			valores_admissiveis_t = []
			lag= lagdistance*(l+1)
			for p in range(0,len(distanciah)):
				if (distanciah[p] < dmaximo):
					limitemin = lag - lineartolerance
					limitemax = lag + lineartolerance
					if (distanciah[p] > limitemin and distanciah[p] < limitemax):					
						if (distanciaxy[p] > 0.000): 				
							check_azimute = (distanciax[p]*cos_Azimute +distanciay[p]*sin_Azimute)/distanciaxy[p]
						else:
							check_azimute = 1
						check_azimute = math.fabs(check_azimute)
						if (check_azimute >= htolerance):							
							check_bandh = (cos_Azimute*distanciay[p]) - (sin_Azimute*distanciax[p])
							check_bandh = math.fabs(check_bandh)
							if (check_bandh < hband):							
								if(distanciah[p] > 0.000):
									check_dip = (math.fabs(distanciaxy[p])*sin_Dip + distanciaz[p]*cos_Dip)/distanciah[p]
								else:
									check_dip = 0.000
								check_dip = math.fabs(check_dip) 
								if (check_dip >= vtolerance):
									check_bandv = sin_Dip*distanciaz[p] - cos_Dip*math.fabs(distanciaxy[p])
									check_bandv = math.fabs(check_bandv)
									if (check_bandv < vband):
										if (check_dip <0 and check_azimute <0):
											valores_admissiveis_h.append(rabo[p])
											valores_admissiveis_t.append(cabeca[p])
										else:
											valores_admissiveis_h.append(cabeca[p])
											valores_admissiveis_t.append(rabo[p])	
			if (len(valores_admissiveis_h) > 0 and len(valores_admissiveis_t) > 0):				
				v_valores_admissiveis_h.append(valores_admissiveis_h)
				v_valores_admissiveis_t.append(valores_admissiveis_t)
		

		for i in range(0,len(v_valores_admissiveis_h)):
			vh = []
			vt = []
			vh = np.array(v_valores_admissiveis_h[i])
			vt = np.array(v_valores_admissiveis_t[i])
			

			slope, intercept, r_value, p_value, std_err = stats.linregress(vh,vt)
			plt.title("Hscatter plot for lag: " + str(lagdistance*(i+1)) + " correlation: " + str(r_value))			
			plt.plot(vh,vt,'ro',vh,slope*vh+intercept,'b-')
			plt.show()
	

		
	
					
    	return True 

    def finalize(self):
        print "Finalize program"
        return True
    def name(self):
        return "hscatter"
##############################################
def get_plugins():
    return["hscatter"]
