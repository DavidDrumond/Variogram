#!/bin/python
import sgems
import matplotlib.pyplot as plt
from matplotlib import colors,ticker,cm
from scipy.interpolate import griddata
from scipy import interpolate
import numpy as np
import math
import random

#..............................................................
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
#..............................................................



"""
VARIOMAP
"""

class   vario_map:
    def _init_(self):
        pass 
    def initialize(self,params):
        self.params = params
        return True
    def execute(self):


	'''
		VARMAP 
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
	C_variogram = int(self.params['C_variogram']['value'])
	C_covariance= int(self.params['C_covariance']['value'])
	nlags= int(self.params['nlags']['value'])
	lagdistance = float(self.params['lagdistance']['value'])
	lineartolerance= float(self.params['lineartolerance']['value'])
	htolerance = float(self.params['htolangular']['value'])
	vtolerance = float(self.params['vtolangular']['value'])
	hband = float(self.params['hband']['value'])
	vband = float(self.params['vband']['value'])
	dip = float(self.params['Dip']['value'])
	azimute = float(self.params['Azimute']['value'])
	gdiscrete = float(self.params['Gdiscrete']['value'])
	ncontour = int(self.params['ncontour']['value'])
	
	# Get values from sgems properties 

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

	# Verify dimensions 
	
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
		


		# Define maximum distance permisible 
		
		dmaximo = (nlags + 0.5)*lagdistance
		
		# Transform tolerances if they are zero 

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
	
		# Define all euclidian distances 
		
		cabeca = []
		rabo = []
		cabeca2 =[]
		rabo2 = []
		distanciah =[]	
		distanciaxy =[]	
		distanciax =[]
		distanciay = []
 		distanciaz = []
		for t in range(0, len(yt)):
			for h in range(t, len(yh)):
				cabeca.append(vh[h])
				rabo.append(vt[t])
				cabeca2.append(vt[h])
				rabo2.append(vh[t])
				dx = xh[h]-xt[t] 				
				dy= yh[h]-yt[t]
				dz = zh[h]-zt[t]
				if (distanciah > 0):
					distanciay.append(dy)
					distanciax.append(dx)
					distanciaz.append(dz)
					distanciah.append(math.sqrt(math.pow(dy,2) + math.pow(dx,2) + math.pow(dz,2)))						
					distanciaxy.append(math.sqrt(math.pow(dy,2) + math.pow(dx,2)))
		
		# Calculate all cosine and sine values 	
		

		cos_Azimute = []
		sin_Azimute = []
		cos_Dip = []
		sin_Dip = []
		flutuante_d = 0
		
		
		for a in range(0,360,10):
			diferenca = a + azimute 
			cos_Azimute.append(math.cos(math.radians(90-diferenca)))
			sin_Azimute.append(math.sin(math.radians(90-diferenca)))
			if (dip != 90 and dip != 180):
				flutuante_d = math.cos(math.radians(90) - math.fabs(math.atan(math.tan(math.radians(dip))*math.sin(math.radians(90-diferenca)))))
				cos_Dip.append(flutuante_d)
				flutuante_d = math.sin(math.radians(90) - math.fabs(math.atan(math.tan(math.radians(dip))*math.sin(math.radians(90-diferenca)))))
				sin_Dip.append(flutuante_d)
			else:	
				cos_Dip.append(1)
				sin_Dip.append(0)
	
		distancias_admissiveis = []
		azimute_admissiveis = []
		dip_admissiveis =[]
		v_valores_admissiveis_h =[]
		v_valores_admissiveis_t =[]
		v_valores_admissiveis_h2 =[]
		v_valores_admissiveis_t2 =[]
		pares = []		

		# Calculate permissible experimental pairs 

		for a in range(0,36):
			azm = math.radians(a*10)
			if (dip != 90 and dip != 180):
				dipatua = math.radians(90) - math.atan(math.tan(math.radians(dip))*math.sin(math.radians(90-diferenca)))
			else:
				dipatua = math.radians(90)
			for l in range(0, nlags):			
				valores_admissiveis_h = []
				valores_admissiveis_t = []
				valores_admissiveis_h2 =[]
				valores_admissiveis_t2 =[]
				distancia_admissivel = []
				azimute_admissivel =[]
				dip_admissivel=[]
				lag= lagdistance*(l+1)
				par= 0 
				for p in range(0,len(distanciah)):
					if (distanciah[p] < dmaximo):
						limitemin = lag - lineartolerance
						limitemax = lag + lineartolerance
						if (distanciah[p] > limitemin and distanciah[p] < limitemax):						
							if (distanciaxy[p] > 0.000): 				
								check_azimute = (distanciax[p]*cos_Azimute[a] +distanciay[p]*sin_Azimute[a])/distanciaxy[p]
							else:
								check_azimute = 1
							check_azimute = math.fabs(check_azimute)
							if (check_azimute >= htolerance):								
								check_bandh = (cos_Azimute[a]*distanciay[p]) - (sin_Azimute[a]*distanciax[p])
								check_bandh = math.fabs(check_bandh)
								if (check_bandh < hband):								
									if(distanciah[p] > 0.000):
										check_dip = (math.fabs(distanciaxy[p])*sin_Dip[a] + distanciaz[p]*cos_Dip[a])/distanciah[p]
									else:
										check_dip = 0.000
									check_dip = math.fabs(check_dip) 
									if (check_dip >= vtolerance):
										check_bandv = sin_Dip[a]*distanciaz[p] - cos_Dip[a]*math.fabs(distanciaxy[p])
										check_bandv = math.fabs(check_bandv)
										if (check_bandv < vband):
											if (check_dip < 0 and check_azimute < 0):
													valores_admissiveis_h.append(rabo[p])
													valores_admissiveis_t.append(cabeca[p])
													valores_admissiveis_h2.append(rabo2[p])
													valores_admissiveis_t2.append(cabeca2[p])
													distancia_admissivel.append(distanciah[p])
													par = par + 1	
													azimute_admissivel.append(azm)
													dip_admissivel.append(dipv)	
											else:
													valores_admissiveis_h.append(cabeca[p])
													valores_admissiveis_t.append(rabo[p])
													valores_admissiveis_h2.append(cabeca2[p])
													valores_admissiveis_t2.append(rabo2[p])
													distancia_admissivel.append(distanciah[p])
													par = par + 1	
													azimute_admissivel.append(azm)
													dip_admissivel.append(dipatua)		
				if (len(valores_admissiveis_h) > 0 and len(valores_admissiveis_t) > 0):	
					if (par > 0):		
						v_valores_admissiveis_h.append(valores_admissiveis_h)
						v_valores_admissiveis_h2.append(valores_admissiveis_h2)
						v_valores_admissiveis_t2.append(valores_admissiveis_t2)
						v_valores_admissiveis_t.append(valores_admissiveis_t)
						distancias_admissiveis.append(distancia_admissivel)
						pares.append(par)
						azimute_admissiveis.append(azimute_admissivel)	
						dip_admissiveis.append(dip_admissivel)
		

		# Calculate continuity functions 	

		# Variogram 
		
		if (C_variogram == 1): 		
			continuidade =[]
			lag_adm =[]
			azimute_adm =[]
			dip_adm =[]
			for i in range(0,len(v_valores_admissiveis_h)):
				flutuantet =[]
				flutuanteh =[]
				flutuanteh2 = []
				flutuantet2 =[]
				flutuanted =[]
				flutuantea=[]
				flutuantedip=[]
				par_var = 0
				flutuanteh = v_valores_admissiveis_h[i][:]
				flutuanteh2= v_valores_admissiveis_h2[i][:]
				flutuantet = v_valores_admissiveis_t[i][:]
				flutuantet2 = v_valores_admissiveis_t2[i][:]
				flutuanted = distancias_admissiveis[i][:]
				flutuantea= azimute_admissiveis[i][:]
				flutuantedip = dip_admissiveis[i][:]
				par_var= pares[i]
				soma = 0
				lagmedio =0
				agmedio =0 
				dgmedio = 0
				for j in range(0, len(flutuanteh)):
					soma = soma + (flutuanteh[j] - flutuantet2[j])*(flutuanteh2[j]-flutuantet[j])/(2*pares[i])
				continuidade.append(soma)
				for z in range(0, len(flutuanted)):
					lagmedio = lagmedio + flutuanted[z]/len(flutuanted)
				lag_adm.append(lagmedio)
				for g in range(0, len(flutuantea)):
					agmedio = agmedio + flutuantea[g]/len(flutuantea)
				azimute_adm.append(agmedio)
				for t in range(0, len(flutuantedip)):
					dgmedio = dgmedio + flutuantedip[t]/len(flutuantedip)
				dip_adm.append(dgmedio)
				
			
			
		# Covariogram 
 
		if (C_covariance == 1): 		
			continuidade =[]
			lag_adm =[]
			azimute_adm =[]
			dip_adm =[]
			for i in range(0,len(v_valores_admissiveis_h)):
				flutuantet =[]
				flutuanteh =[]
				flutuanted =[]
				flutuantea=[]
				flutuantedip=[]
				par_var = 0
				flutuanteh = v_valores_admissiveis_h[i]
				flutuantet = v_valores_admissiveis_t[i]
				flutuanted = distancias_admissiveis[i]
				flutuantea= azimute_admissiveis[i]
				flutuantedip = dip_admissiveis[i]
				par_var= pares[i]
				soma = 0
				lagmedio =0
				agmedio =0 
				dgmedio = 0
				somah = 0
				somat = 0 
				mediah = 0
				mediat =0
				for d in range (0, len(flutuanteh)):
					somah = somah + flutuanteh[d]
				mediah = float(somah/len(flutuanteh))
				for t in range (0, len(flutuantet)):
					somat = somat + flutuantet[t]
				mediat = float(somat/len(flutuantet))
				for j in range(0, len(flutuanteh)):		
						soma = soma + float(((flutuanteh[j] - mediah)*(flutuantet[j] - mediat))/(par_var))
				continuidade.append(soma)
				for z in range(0, len(flutuanted)):
					lagmedio = lagmedio + flutuanted[z]/len(flutuanted)
				lag_adm.append(lagmedio)
				for g in range(0, len(flutuantea)):
					agmedio = agmedio + flutuantea[g]/len(flutuantea)
				azimute_adm.append(agmedio)
				for y in range(0, len(flutuantedip)):
					dgmedio = dgmedio + flutuantedip[y]/len(flutuantedip)
				dip_adm.append(dgmedio)
				
		
		
		
		
		# Plot variograms on map 
	
		
		x = np.array(lag_adm)*np.sin(np.array(azimute_adm))
		y = np.array(lag_adm)*np.cos(np.array(azimute_adm))



		Xi = np.linspace(-max(x)-lagdistance,max(x)+lagdistance,gdiscrete)
		Yi = np.linspace(-max(y)-lagdistance,max(y)+lagdistance,gdiscrete)

	

		#make the axes
		f = plt.figure()
		left, bottom, width, height= [0,0.1, 0.7, 0.7]
		ax  = plt.axes([left, bottom, width, height])
		pax = plt.axes([left, bottom, width, height],
				projection='polar',
				axisbg='none')
		
		pax.set_theta_zero_location("N")
		pax.set_theta_direction(-1)

		cax = plt.axes([0.8, 0, 0.05, 1])
		ax.set_aspect(1)
		ax.axis('Off')
	

		# grid the data.
		Vi = griddata((x, y), np.array(continuidade), (Xi[None,:], Yi[:,None]), method='cubic')	
		cf = ax.contourf(Xi,Yi,Vi, ncontour, cmap=plt.cm.jet)


		gradient = np.linspace(-max(continuidade),max(continuidade), math.fabs(2*max(continuidade)))
		gradient = np.vstack((gradient, gradient))
		cax.xaxis.set_major_locator(plt.NullLocator())
		cax.yaxis.tick_right()
		cax.imshow(gradient.T, aspect='auto', cmap=plt.cm.jet)

		plt.show()
		plt.close()
		
	
    def finalize(self):
        print "Finalize program"
        return True
    def name(self):
        return "vario_map"
##############################################
def get_plugins():
    return["vario_map"]
