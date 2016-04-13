#!/bin/python
import sgems
import matplotlib.pyplot as plt
from matplotlib import colors,ticker,cm
from scipy.interpolate import griddata
from scipy import interpolate
import numpy as np
import math
import random



"""
VARIOMAP
"""

#.................................................................
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
#....................................................................		




class   variogram:
    def _init_(self):
        pass 
    def initialize(self,params):
        self.params = params
        return True
    def execute(self):
	print(self.params)

	'''
		EXPERIMENTAL VARIOGRAMS CALCULATION 
		_____________________________________________________________________________

		This is a template for SGems program for calculating experimental variograms.
		The user have to fill the parameters with a ui.file with same name of this python 
		file before start program. The QT userinterface connect SGems data with this routine
		using the class params.
		
		The output of this program is two files, one relatory with all experimental points 
		used in Automatic fitting procedure and other with a hml file used to import 
		experimental variograms in SGems

	AUTHOR: DAVID ALVARENGA DRUMOND 		2016
	
	'''	

	'''
	INITIAL PARAMETERS 
	______________________________________________________________________________

	All initial parameters are transcribe by ui.file interface. Variables are choose 
	among common variables. The follow variables are:

	COMMON VARIABLES: 

	head_property= adress of head property
	head_grid = adress of head grid 
	tail_property = adress of tail property
	tail_grid= adress of tail grid 
	C_variogram = check for variogram function calculation 
	C_covariance= check for covariance function calculation 
	C_relative_variogram = check for relative variogram calculation
	C_madogram = check for madogram calculation 
	C_correlogram = check for correlogram calculation 
	C_PairWise = check for pairwise calcluation 
	nlags= number of lags in experimental variogram 
	lagdistance = lenght of lag in experimental variogram 
	lineartolerance= linear tolerance of experimental variogram 
	htolerance = horizontal angular tolerance
	vtolerance = vertical angular tolerance
	hband = horizontal bandwidth
	vband = vertical bandwidth
	nAzimuths = number of azimuths 
	NDips = number of dips 
	Azimth_diference = angular diference between azimuths 
	Dip_difference = angular diference between dips
	sAzimuth = initial azimuth value
	sDip = initial dip value 
	inv = choose of invert correlogram axis 
	save = save adress of relatory file 
	save2 = save adress of Sgems variogram file 
	nhead = number of head to print in relatory file 
	ntail = number of tail to print in relatory file 

	
	'''

	# Stabilish initial parameters 

	head_property= self.params['head_prop']['property']
	head_grid = self.params['head_prop']['grid']
	tail_property = self.params['tail_prop']['property']
	tail_grid= self.params['tail_prop']['grid']
	C_variogram = int(self.params['C_variogram']['value'])
	C_covariance= int(self.params['C_covariance']['value'])
	C_relative_variogram = int(self.params['C_relative_variogram']['value'])
	C_madogram = int(self.params['C_madogram']['value'])
	C_correlogram = int(self.params['C_correlogram']['value'])
	C_PairWise = int(self.params['C_PairWise']['value'])
	nlags= int(self.params['nlags']['value'])
	lagdistance = float(self.params['lagdistance']['value'])
	lineartolerance= float(self.params['lineartolerance']['value'])
	htolerance = float(self.params['htolangular']['value'])
	vtolerance = float(self.params['vtolangular']['value'])
	hband = float(self.params['hband']['value'])
	vband = float(self.params['vband']['value'])
	nAzimuths = int(self.params['nAzimuths']['value'])
	NDips = int(self.params['NDips']['value'])
	Azimth_diference = float(self.params['Azimth_diference']['value'])
	Dip_difference = float(self.params['Dip_difference']['value'])
	sAzimuth = float(self.params['sAzimuth']['value'])
	sDip = float(self.params['sDip']['value'])
	inv = int(self.params['Inv']['value'])
	save = self.params['Save']['value']
	save2 = self.params['Save2']['value']
	nhead = self.params['NHEAD']['value']
	ntail = self.params['NTAIL']['value']
	
	# Obtain properties in SGems
	

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

	# Verify number of dimensions in the problem 
	
	# Verify number of dimensions in head 

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

	# Verify number of dimensions in tail 

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
		


		# Define maximum permissible distance
		
		dmaximo = (nlags + 0.5)*lagdistance
		
		#Transform tolerances if they are zero 

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
		
		
		#Convert tolerances in radians 
		htolerance = math.radians(htolerance)
		vtolerance = math.radians(vtolerance)	

		# Determine projections of Dip and Azimuth 
		htolerance = math.cos(htolerance)
		vtolerance = math.cos(vtolerance)
	
		# Define all direction distances 
		
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

		# Calculate all sin and cosine functions of Dip and Azimuth 	
		

		
		azimute = []
		fAzimuth = float(Azimth_diference*(nAzimuths-1)) + sAzimuth 
		azimute = np.linspace(sAzimuth,fAzimuth,nAzimuths)
		

		dip = [] 
		fDip = Dip_difference*(NDips-1) + sDip
		dip = np.linspace(sDip,fDip,NDips)
	

		cos_Azimute=[]
		sin_Azimute=[]
		for a in azimute:
			cos_Azimute.append(math.cos(math.radians(90-a)))
			sin_Azimute.append(math.sin(math.radians(90-a)))
		cos_Dip =[]
		sin_Dip =[]
		for a in dip:	
			cos_Dip.append(math.cos(math.radians(90-a)))
			sin_Dip.append(math.sin(math.radians(90-a)))
	
		
		

		distancias_admissiveis = []
		azimute_admissiveis = []
		dip_admissiveis =[]
		v_valores_admissiveis_h =[]
		v_valores_admissiveis_t =[]
		v_valores_admissiveis_h2 =[]
		v_valores_admissiveis_t2 =[]
		pares = []		

		# Calculate admissible pairs  
		for a in range(0,len(azimute)):
			azm = azimute[a]
			for d in range(0, len(dip)):
				dipv = dip[d]
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
											check_dip = (math.fabs(distanciaxy[p])*sin_Dip[d] + distanciaz[p]*cos_Dip[d])/distanciah[p]
										else:
											check_dip = 0.000
										check_dip = math.fabs(check_dip) 
										if (check_dip >= vtolerance):
											check_bandv = sin_Dip[d]*distanciaz[p] - cos_Dip[d]*math.fabs(distanciaxy[p])
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
													dip_admissivel.append(dipv)
					if (len(valores_admissiveis_h) > 0 and len(valores_admissiveis_t) > 0):			
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
				
		# Correlogram
		if (C_correlogram == 1): 		
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
				diferencah =0
				diferencat =0
				mediah = 0
				mediat =0
				varh = 0
				vart = 0
				for d in range (0, len(flutuanteh)):
					somah = somah + flutuanteh[d]
				mediah = float(somah/len(flutuanteh))
				for t in range (0, len(flutuantet)):
					somat = somat + flutuantet[t]
				mediat = float(somat/len(flutuantet))
				for u in range(0, len(flutuanteh)):
					diferencah = flutuanteh[u]**2 + diferencah
				varh = math.sqrt(float(diferencah)/len(flutuanteh)-mediah**2)
				for m in range(0, len(flutuanteh)):
					diferencat = flutuantet[m]**2 + diferencat				
				vart = math.sqrt(float(diferencat)/len(flutuantet)-mediat**2)
				for j in range(0, len(flutuanteh)):
						if (vart > 0 and varh >0 ):
							if (inv == 1):
								soma = soma + 1-(flutuanteh[j] - mediah)*(flutuantet[j] - mediat)/(par_var*vart*varh)
							else:
								soma = soma + (flutuanteh[j] - mediah)*(flutuantet[j] - mediat)/(par_var*vart*varh)	
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

		# PairWise 
		if (C_PairWise == 1): 		
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
				flutuanteh = v_valores_admissiveis_h[i][:]
				flutuantet = v_valores_admissiveis_t[i][:]
				flutuanted = distancias_admissiveis[i][:]
				flutuantea= azimute_admissiveis[i][:]
				flutuantedip = dip_admissiveis[i][:]
				par_var= pares[i]
				soma = 0
				lagmedio =0
				agmedio =0 
				dgmedio = 0
				for j in range(0, len(flutuanteh)):
					s = ((flutuanteh[j] - flutuantet[j])/(flutuanteh[j]+flutuantet[j]))**2
					soma = soma + 2*s/(1.0*par_var)
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
				
	
		# Relative Variogram
		if (C_relative_variogram  == 1): 		
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
				flutuanteh = v_valores_admissiveis_h[i]
				flutuantet = v_valores_admissiveis_t[i]
				flutuanteh2 = v_valores_admissiveis_h2[i]
				flutuantet2 = v_valores_admissiveis_t2[i]
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
						soma = soma + float((flutuanteh[j] -flutuantet2[j])*(flutuanteh2[j] - flutuantet2[j])/(par_var*(mediat+mediah)**2/4))
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
		
		# Madogram 
		if (C_madogram == 1): 		
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
				for j in range(0, len(flutuanteh)):
						soma = soma + float(math.fabs(flutuanteh[j] - flutuantet[j])/(2*par_var))
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
		

		# Write experimental variograms in relatory 


		save= open(save,"a")
		save.write(nhead + "\n")
		save.write(ntail + "\n")
		save.write(" Azimuth	Dip	lag	variogram	pairs \n")
		for i in range(0, len(continuidade)):
			save.write(str(azimute_adm[i]) + "	" +str(dip_adm[i]) + "	" + str(lag_adm[i]) + "	" + str(continuidade[i])+ "	" + str(pares[i]) + "\n") 	

		# Separate experimental variograms according prescribe direction 

		posicao =[]
		posicao.append(-1)
		for i in range(1, len(azimute_adm)):
			if (round(azimute_adm[i],2) != round(azimute_adm[i-1],2) or round(dip_adm[i],2) != round(dip_adm[i-1],2)):
				posicao.append(i-1)	
		posicao.append(len(azimute_adm)-1)
		
		azimutes = []
		dips = []
		lags =[]
		continuidades =[]
		paresv =[]

		for i in range(1,len(posicao)):
			azimutes.append(azimute_adm[posicao[i-1]+1:posicao[i]])
			dips.append(dip_adm[posicao[i-1]+1:posicao[i]])
			lags.append(lag_adm[posicao[i-1]+1:posicao[i]])
			continuidades.append(continuidade[posicao[i-1]+1:posicao[i]])
			paresv.append(pares[posicao[i-1]+1:posicao[i]])

		# Write hml file for experimental variograms

		save2= open(save2,"w")
		save2.write("    <experimental_variograms> \n")

		for i in range(0, len(azimutes)):
			azim = []
			Dip =[]
			lagv =[]
			cont =[]
			par =[]
			azim = azimutes[i]
			Dip = dips[i]
			lagv = lags[i]
			cont =continuidades[i]
			par = paresv[i]
			save2.write("  <variogram> \n")
			save2.write("  <title>variogram -azth="+str(round(azim[0],2))+", dip="+str(round(Dip[0],2))+"</title> \n")
			direction1 = math.sin(Dip[0])*math.cos(azim[0])
			direction2 = math.sin(Dip[0])*math.sin(azim[0])
			direction3 = -math.sin(Dip[0])
			save2.write("<direction>" +str(direction1) + " " + str(direction2) +" "+ str(direction3)+ "</direction> \n")
			save2.write("   <x>")
			for j in range(0,len(lagv)):
				save2.write(str(lagv[j]) + " ")
			save2.write("</x> \n")
			save2.write("   <y>")
			for p in range(0,len(cont)):
				save2.write(str(cont[p])+ " ")
			save2.write("</y> \n")
			save2.write("   <pairs>")
			for h in range(0,len(cont)):
				save2.write(str(par[h])+ " ")
			save2.write("</pairs> \n")
			save2.write("  </variogram> \n")
		save2.write("   </experimental_variograms> \n")
		save2.close()	
			
				
	
	print("Finish")
			
    	return True 

    def finalize(self):
        print "Finalize program"
        return True
    def name(self):
        return "variogram"
##############################################
def get_plugins():
    return["variogram"]
