#!/bin/python
import sgems
import matplotlib.pyplot as plt
import numpy as np
import math
import random



"""
2016 - PORTO ALEGRE - UNIVERSIDADE FEDERAL DO RIO GRANDE DO SUL 
PROGRAMA FREEWARE PARA DESENVOLVIMENTO DO SOFTWARE SGEMS
PROGRAMA PARA MODELAGEM DE VARIOGRAMAS AUTOMATICOS EM UM MODELO LINEAR DE CORREGIONALIZACAO 

"""


#...............................................................................

# Function for calculating weights of variogram number of pairs  
def pesos_pares(par):
    soma =0
    pesos = []
    for i in range(0,len(par)):
        soma = soma + par[i]
    for j in range(0, len(par)):
        pesos.append(float(par[j])/soma)
    return pesos

# ..............................................................................

# Function for calculating weights of variogram lag distance
def pesos_distancia(distancia):
    soma = 0
    pesos =[]
    for i in distancia:
        inverso = 1/i
        soma = soma + inverso
    for i in distancia:
        pesos.append((1/i)/soma)
    return pesos

#...............................................................................

# Fuction to calculate variogram model values for specify parameters
def modelo_variogama(range_estruturas, sill_estruturas, modelos_estruturas, lags, tamanho_lags, tamanho_estruturas):

    valor_variograma_estrutura = []
    n_estruturas = len(modelos_estruturas)
    n_lags = len(lags)
    for i in range(0, n_estruturas):
        valor_variograma = []
        for j in range(0, n_lags):
            if modelos_estruturas[i] == 0:
                if lags[j] < range_estruturas[i]:
                    valor_variograma.append(sill_estruturas[i]*(1.5*lags[j]-0.5*pow((lags[j]/range_estruturas[i]),3)))
                else:
                    valor_variograma.append(sill_estruturas[i])
            elif modelos_estruturas[i] == 1:
                valor_variograma.append(sill_estruturas[i]*(1-math.exp(-3*lags[j]/range_estruturas[i])))
            elif modelos_estruturas[i] == 2:
                valor_variograma.append(sill_estruturas[i]*(1-math.exp(-3*pow((lags[j]/range_estruturas[i]),2))))
        valor_variograma_estrutura.append(valor_variograma)
    return valor_variograma_estrutura

#...............................................................................

# Function to calculate variogram modelling residuals 
def desvio_quad(variograma, modelo_variograma, tamanho_variograma, tamanho_estruturas, efeito_pepita, pesos_distancia, pesos_pares):
    soma_desvio = 0
    for i in range(0, tamanho_estruturas):
        diferenca = 0
        m = list(modelo_variograma[i])
        for j in range(0, tamanho_variograma):
            diferenca = pesos_distancia[j]*pesos_pares[j]*math.pow((m[j]-variograma[j]-efeito_pepita),2)/tamanho_variograma + diferenca
        soma_desvio = soma_desvio + diferenca
        soma_desvio = math.fabs(soma_desvio)
    return soma_desvio
#....................................................................................



# ------------------------Start plugin class ----------------------------------------------------- 


class   automatic_fitting2:
    def _init_(self):
        caminho = "DADOS"
        numero_variogramas = 1
        numero_estruturas = 1
        azimute_direcao_principal = 0
        dip_direcao_principal = 0
        numero_interacoes = 0
        pass 
    def initialize(self,params):
        self.params = params
        return True
    def execute(self):

	#----------------------------------------------------------------------------------------

	'''
		AUTOMATIC VARIOGRAM MODELLING USING A LINEAR MODEL OF CORREGIONALIZATION
		_____________________________________________________________________________

		This is a template for SGems program for calculating automatic variograms using 
		a linear model of corregionalization. The user have to fill the parameters with 
		a ui.file with same name of this python file before start program. The QT user 
		interface connect SGems data with this routine using the class params.
		
		The output of this program is a graphic with all variogram models fitted and 
		lines printed in command prompt with variogram parameters

	AUTHOR: DAVID ALVARENGA DRUMOND 		2016
	
	'''	

	'''
	INITIAL PARAMETERS 
	______________________________________________________________________________

	All initial parameters are transcribe by ui.file interface. Variables are choose 
	among common variables in automatic fitting and restriction variables used for 
	marginal conditions. The follow variables are:

	COMMON VARIABLES: 

	caminho = File adress of experimental variogram files 
	numero_variogramas = number of variograms 
	numero_estruturas = number of structures 
	numero_interacoes = number of interations of monte carlo algorithm
	set_maximumrange = variable for choose maximum range optimization 
	set_minrange = variable for choose minimum range optimization 
	set_verticalrange = variable for choose vertical range optimization 
	azimute_direcao_principal = principal direction azimuth 
	dip_direcao_principal = principal direction dip 
	azimute_direcao_secundaria = secondary direction azimuth 
	dip_direcao_secundaria = secondary direction azimuth 
	azimute_direcao_vertical = vertical direction azimuth 
	dip_direcao_vertical = vertical direction dip 
	restriction = choose of restrictions in modelling 
	nlags = number of lags in each variogram 
	nvariables = number of variables used in modelling 
	min_npairs = minimum number of pairs in variogram experimental to use in modelling

	RESTRICTION VARIABLES:

	min_contribution = minimum variogram contribution 
	max_contribution = maximum variogram contribution 
	min_nugget = minimum nugget effect 
	max_nugget = maximum nugget effect 
	min_range_max = minimum range of direction of maximum continuity
	max_range_max = maximum range of direction of maximum continuity 
	min_range_min = minimum range of direction of secondary continuity
	max_range_min = maximum range of direction of secondary continuity 
	min_range_vert = minimum range of direction of vertical continuity 
	max_range_vert = maximum range of direction of vertical continuity 
	
	'''
        caminho = str(self.params['endereco_arquivo']['value'])
        numero_variogramas = int(self.params['N_variogramas']['value'])
        numero_estruturas = int(self.params['N_estruturas']['value'])
        numero_interacoes =int(self.params['N_interacoes']['value'])
	
	set_maximumrange = int(self.params['set_maximumrange']['value'])
	set_minrange = int(self.params['set_minrange']['value'])
	set_verticalrange = int(self.params['set_verticalrange']['value'])
	
        if (set_maximumrange ==1):
		azimute_direcao_principal =float(self.params['Azimute_dp']['value'])
        	dip_direcao_principal = float(self.params['Dip_dp']['value'])

	if (set_minrange == 1):
		azimute_direcao_secundaria = float(self.params['Azimute_dm']['value'])
		dip_direcao_secundaria = float(self.params['Dip_dm']['value'])
	
	if (set_verticalrange == 1):
		azimute_direcao_vertical = float(self.params['Azimute_dv']['value'])
		dip_direcao_vertical = float(self.params['Dip_dv']['value'])

	restriction = int(self.params['restriction']['value'])
	nlags = int(self.params['nlags']['value'])
	nvariables = int(self.params['nvariables']['value'])
	min_npairs = int(self.params['min_npairs']['value'])
	
	if (restriction == 1):	
		min_contribution = float(self.params['min_contribution']['value'])
		max_contribution = float(self.params['max_contribution']['value'])
		min_nugget = float(self.params['min_nugget']['value'])
		max_nugget = float(self.params['max_nugget']['value'])
		if (set_maximumrange == 1):
			min_range_max = float(self.params['min_range_max']['value'])
			max_range_max = float(self.params['max_range_max']['value'])
		if (set_minrange ==1):
			min_range_min = float(self.params['min_range_min']['value'])
			max_range_min = float(self.params['max_range_min']['value'])
		if (set_verticalrange == 1):
			min_range_vert = float(self.params['min_range_vert']['value'])
			max_range_vert = float(self.params['max_range_vert']['value'])
		
	else:
		min_contribution = 0
		max_contribution = 0
		min_nugget = 0
		max_nugget =0
		min_range_max = 0
		max_range_max = 0
		min_range_min = 0
		max_range_min = 0
		min_range_vert = 0
		max_range_vert = 0

	# ----------------------------------------------------------------------------------------

	'''
	ROUTINE FOR FILE READING 
	__________________________________________________________________________________
	
	This routine is used for reading the experimental variogram file content. 	

	1) Test if experimental variogram file can be open 
	2) Read variogram file and put experimental variogram in lists 
		a) azimuteVV = azimuth of experimental variograms 
		b) dipVV = dip of experimental variograms 
		c) variogramasV = experimental variogram values 
		d) lagsV = lags of each experimental variograms 
		e) paresV = number of pairs in each experimental variograms 
		f) indice_head = index of head values variable 
		g) indice_tail = index of tail values variable
	3) Determine index for variograms of principal direction, secondary and vertical. Filter 
	which experimental variograms will be use. 

	''' 


        # 1) Test Variogram file open 
        try:
		arquivo = open (caminho, "r")
        except:
                print("Could not open the file! Please try again!")
        
	size_of_matriz = math.sqrt(nvariables)

	# 2) Read Variogram files 	

	azimuteVV =[]
	dipVV =[]
	variogramasV =[]
	dipsV =[]
	lagsV =[]
	paresV=[]
	indice_head =[]
	indice_tail =[]
	tipo = []
	
	
	for p in range(0,nvariables):	
			
		# Read index of head values of variogram 
		linha = arquivo.readline() 
		indice_head.append(int(linha))       
		
		# Read index of tail values of variogram
	        linha = arquivo.readline()
		indice_tail.append(int(linha))
		
		# Read header file 
		linha = arquivo.readline()

	        # Define vectors for variogram parameters
	        azimuteV = []
	        dipV = []
	        variogramas =[]
	        dips = []
	        lags = []
	        pares =[]

	        # For each variogram do 
	        for i in range(0,numero_variogramas):
		    variograma =[]
	            lag = []
	            par =[]
		    for j in range(0, nlags):
	            	linha = arquivo.readline()
	            	linhasplit = linha.split()
			if(int(linhasplit[4]) >= min_npairs):
		    		if (j ==0):
					azimuteV.append(float(linhasplit[0]))
					dipV.append(float(linhasplit[1]))
		    		lag.append(float(linhasplit[2]))
		    		variograma.append(float(linhasplit[3]))
		    		par.append(int(linhasplit[4]))
		    variogramas.append(variograma)
		    lags.append(lag)
		    pares.append(par)
		azimuteVV.append(azimuteV)
		dipVV.append(dipV)
		variogramasV.append(variogramas)
		lagsV.append(lags)
		paresV.append(pares)

        arquivo.close()

        # 3) Determine variogram index for principal, secondary and vertical directions

        I_dir_principal = 0
	I_dir_secundario = 0
	I_dir_vertical = 0

	informado_principal = False
	informado_secundario = False
	informado_vertical = False

        indice = 0
	interacoes = len(azimuteV)
        for i in range(0,interacoes):
	    if (set_maximumrange ==1):
		    if (azimuteV[i] == azimute_direcao_principal and dipV[i] == dip_direcao_principal):
		    	I_dir_principal = indice
			informado_principal = True
	    if (set_minrange == 1):
		    if (azimuteV[i] == azimute_direcao_secundaria and dipV[i] == dip_direcao_secundaria):
		    	I_dir_secundario = indice
			informado_secundario = True
	    if (set_verticalrange == 1):  
		    if (azimuteV[i] == azimute_direcao_vertical and dipV[i] == dip_direcao_vertical):
		    	I_dir_vertical = indice
			informado_vertical = True
	    indice = indice + 1

	# ------------------------------------------------------------------------------------------	
	
	'''

	START FUNCTION OF VARIOGRAM MODELLING 
	_______________________________________________________________________________

	This fucntion optimizing variogram modelling for each direction specify in last procedure.
	It optimizes variogram principal direction and copy some parameters of principal direction 
	in secondary and vertical directions. 
	this fucntion requires as input the following sequence:

	a) I_dir_principal = index of direction to model 
	b) variogramasV = experimental variograms in this direction 
	c) lasgV = experimental variogram distances in this direction 
	d) paresV = experimental variogram pairs in this direction 
	e) nvariables = number of variables in experimental variograms 
	f) numero_estruturas = number of variogram structures previous defined 
	g) restriction = Choose of resctritions in variogram modelling 
	h) min_contribution, max_contribution, min_nugget, max_nugget, min_range, max_range=
		restrictions in variogram modelling 
	i) modeling = Variogram model of principal direction 	
	j) direcao = Direction which is modelled ( 1 principal, 2 secondary, 3 vertical )
	k) pepm = nugget effect of principal direction 
	l) sillm = contribution of principal direction 
	
	The algorithm following the sequence of steps:
	
	1) Defining weights for experimental variogram pairs 
	2) Determining first variogram parameters 
	3) Define sorted parameters 
	4) Sort parameters to model 
 	5) Accept or reject modification based on lower residual
		
	'''		
	def otimizacao(I_dir_principal, variogramasV, lagsV, paresV, nvariables, 		          numero_estruturas,restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range,max_range, modeling, direcao, pepm, sillm): 
		count = 0
		check_LMC = False
		check_LMC2 = False
		while (check_LMC == False and check_LMC2 == False):
			check_LMC = False
			check_LMC2 = False
			residuoV =[]
			rangeV =[]
			sillV =[]
			testerV =[]
			pepV = [] 

			for n in range(0,nvariables):

				variogramas = variogramasV[n]
				pares = paresV[n]
				lags = lagsV[n]
		
				# 1) Define experimental variogram weigths		
	
				p_pares = []
				p_pares = pesos_pares(pares[I_dir_principal][:])
				p_distancia =[]
				p_distancia = pesos_distancia(lags[I_dir_principal][:])
			
				# 2) Determine first variogram model 
				#..................................
				# a) Define sill for each structure
				sill_estruturas = []       
				sill_medio = 0
				soma = 0
				n = 0
				resultado_pepita = 0         
				for i in variogramas[I_dir_principal]:
				    soma = soma + i
				    n = n + 1
				sill_medio = float(soma/(n*numero_estruturas))
				if (count > 1):
					sill_medio = float(soma/(n*numero_estruturas))
				for i in range(0,numero_estruturas):
				    sill_estruturas.append(sill_medio )				
				#b) Define nugget effect 
				resultado_pepita = 0.1*sill_medio
				#c) Define range of each strucutre
				range_estruturas = []
				max_range_2 = float(max(lags[I_dir_principal])/2)
				for i in range(0, numero_estruturas):
				    range_estruturas.append(float(max_range_2/numero_estruturas))
				#d) Define variogram model for each structure
				tipo_estrutura = []
				for i in range(0,numero_estruturas):
				    tipo_estrutura.append(0)
				n = 0
	
				#....................................................................
			
				# 2) Define variogram parameters for direction prescribed 

				variograma = []	       
				variograma = variogramas[I_dir_principal]
				lag = []     
				lag = lags[I_dir_principal]				
				numero_variogramas = len(variograma)
				numero_lags = len(lag)

				#....................................................................

				# 3) Define the number of parameters to random sort 		
				numero_de_parametros = numero_estruturas*4

				
				# a) Define first model based on experimental values 
	
				modelo_estruturas = []	
				modelo_estruturas = modelo_variogama(range_estruturas[:], sill_estruturas[:], tipo_estrutura[:], lag[:], numero_lags, numero_estruturas)    			
				# b) Calculate first residual 
		
				residuo = 0	
				residuo = desvio_quad(variograma, modelo_estruturas[:], numero_lags, numero_estruturas, resultado_pepita,p_distancia[:],p_pares[:])

				# c) Define list of modified parameters

				range_estruturas_modficado = []
				sill_estruturas_modficado = []
				tipo_estruturas_modificado =[]
				pepita_modificado = 0
	

				# d) Copy intial values to modified parameters
		 
				range_estruturas_modificado = range_estruturas[:]
				sill_estruturas_modificado = sill_estruturas[:]
				tipo_estruturas_modificado = tipo_estrutura[:]
				pepita_modificado = resultado_pepita


				# e) Define output variables 

				resultado_residuo = 0
				resultado_range = []
				resultado_estruturas =[]
				resultado_sill = []


				# f) Define convergence of solution 
				convergencia_solucao = False

				#  4) Sorted parameters to model  
				#........................................................................................
				
				for i in range(0,numero_interacoes):

				    # a) Define a random sort parameter to modify 
				    aleatorio_parametro = random.randint(0,numero_de_parametros)

				    # b) Define if increment is positive or negative 
				    aleatorio_sinal  = random.choice([-1,1])

				    # c) According parameter sorted modify itself to 7,5%

				    if (aleatorio_parametro < numero_estruturas):
					range_estrutura = range_estruturas[aleatorio_parametro]
					range_estrutura = range_estrutura + 0.075*aleatorio_sinal*range_estrutura
					if(restriction == 1):
						if (range_estrutura < min_range):				
							range_estrutura = min_range
						if (range_estrutura > max_range):
							range_estrutura = max_range
					range_estruturas_modificado[aleatorio_parametro] = range_estrutura
				    if (aleatorio_parametro >= numero_estruturas and aleatorio_parametro < 2*numero_estruturas):
					sill_estrutura = sill_estruturas[(aleatorio_parametro-numero_estruturas)]
					sill_estrutura = sill_estrutura + 0.075*aleatorio_sinal*sill_estrutura
					if (restriction == 1):
						if(sill_estrutura > max_contribution):
							sill_estrutura = max_contribution
						if(sill_estrutura < min_contribution):
							sill_estrutura = min_contribution
					sill_estruturas_modificado [(aleatorio_parametro - numero_estruturas)] = sill_estrutura

				    # d) If model is sorted choose variogram model despite 0 -Spherical , 1- Exponential, 2-Gaussian
				    if (aleatorio_parametro >= 2*numero_estruturas and aleatorio_parametro < 3*numero_estruturas):
					aleatorio_modelo = random.randint(0,2)
					tipo_estruturas_modificado[aleatorio_parametro - 2*numero_estruturas] = aleatorio_modelo
					
				    # e) Define the total contribution
				    sill_total = 0
				    for j in range(0, numero_estruturas):
					sill_total = sill_total + sill_estruturas_modificado[j]

				    # f) Define a nugget effect if this parameter is choiced 
				    if (aleatorio_parametro >= 3*numero_estruturas and aleatorio_parametro < 4*numero_estruturas):
					pepita_modificado = math.fabs(pepita_modificado + 0.075*aleatorio_sinal*pepita_modificado)
					if (restriction == 1):
						if (pepita_modificado < min_nugget):
							pepita_modificado = min_nugget
						if (pepita_modificado > max_nugget):
							pepita_modificado = max_nugget
									   
				    # g) Calculate the total sill 
		    
			 	    sill_total = sill_total + pepita_modificado
		    
				    range_estruturas_modificado.sort()
			   	    sill_estruturas_modficado.sort()

				    # 5) Define variogram model and choose to accept or reject modification 

				    # a) calculate variogram model 
				    modelo_estruturas = modelo_variogama(range_estruturas_modificado[:], sill_estruturas_modificado[:], tipo_estruturas_modificado[:], lag[:], numero_lags, numero_estruturas)
				    
  				    # b) calculate variogram residual 
				    residuo_modificado = desvio_quad(variograma,modelo_estruturas, numero_lags, numero_estruturas, pepita_modificado, p_distancia, p_pares)

				    # c) Accept variogram parameters modification based on lower residual 
				    if (residuo_modificado < residuo):
					convergencia_solucao = True
					r_range = []
					r_sill = []
					r_testr = []
					r_residuo = residuo_modificado
					r_range = range_estruturas_modificado[:]
					r_sill = sill_estruturas_modificado[:]
					r_testr = tipo_estruturas_modificado[:]		                
					r_pep = pepita_modificado
					residuo = residuo_modificado
					
				if (convergencia_solucao == True):
				    residuoV.append(r_residuo)
				    rangeV.append(r_range)
				    sillV.append(r_sill)
				    testerV.append(r_testr)
				    pepV.append(r_pep)

			
			#........................................................................................
				
			
			'''
			Test LMC model 
			____________________________________________________________________________

			This routine test LMC model if is acceptable. The contribution matrixes are 
			created together with nugget matrixes. Cross variograms are normalizated to 
			express the same model function. The following steps are

			1) Create a contribution matrix 
			2) Regularize contribution matrix with cross variograms 
			3) Calculate determinants for each variogram structure 
			4) Verify if all determinants are greater than zero
			5) Create matrix of nugget effects 
			6) Regularize nugget effect matrix with cross variograms 
			7) Calculate determinants of nugget effect and check if LMC is found	

			'''
			
			#1) Create contribution matrix 
			indice = int(math.sqrt(nvariables))		
			matriz_contr = []
			determ = []
			sillV_new = []
			for p in range(0,numero_estruturas):
				matriz = np.zeros((indice,indice))			
				for x in range(0,len(sillV)):
					
					sill_t = []
					sill_v = []				
					sill_t = sillV[x]
					sill_v = sill_t[p]
				
					i = indice_head[x]						
					j = indice_tail[x]
					i = i - 1
					j = j - 1					
					
					matriz[i][j] = sill_v
				
				# 2) Regularize contribution matrix with cross variograms

				matriz_new = np.zeros((indice,indice))
				for i in range(0,int(math.sqrt(nvariables))):
					for j in range(0,int(math.sqrt(nvariables))):
						if (i==j):
							matriz_new[i][j] = matriz[i][j]
						else:
							matriz_new[i][j] = (matriz[i][j] + matriz[j][i])/2
				sill_t_new =[]
				for i in range(0,int(math.sqrt(nvariables))):
						for j in range(0,int(math.sqrt(nvariables))):
							sill_t_new.append(matriz_new[i][j])
				sillV_new.append(sill_t_new)

				# 3) Calculate determinants for each variogram structure 
		
				determinante = np.linalg.det(matriz_new)
				determ.append(determinante)
				matriz_contr.append(matriz_new)

	
			for p in range(0,numero_estruturas):
				for x in range(0,nvariables):
					sillV[x][p] = sillV_new[p][x]
		
			# 4) Verify if all determinants are greater than zero 

			def verif(v):
				for i in v:
					if i>0:
						pass
					else:
						return False
				return True
			
			if (determ > 0):
				check_LMC = verif(determ)
			

			# 5) Create matrix of nugget effects 
				
			matriz2 = np.zeros((indice,indice))			
			for x in range(0,len(pepV)):
				
				pep_t = []
				pep_v = []				
				pep_v = pepV[x]
				
			
				i = indice_head[x]						
				j = indice_tail[x]
				i = i - 1
				j = j - 1					
					
				matriz2[i][j] = pep_v

			# 6) Regularize nugget effect matrix with cross variograms 
			
			matriz2_new = np.zeros((indice,indice))
			for i in range(0,int(math.sqrt(nvariables))):
				for j in range(0,int(math.sqrt(nvariables))):
					if (i==j):
						matriz2_new[i][j] = matriz2[i][j]
					else:
						matriz2_new[i][j] = (matriz2[i][j] + matriz2[j][i])/2			
			pep_v_new =[]
			for i in range(0,int(math.sqrt(nvariables))):
				for j in range(0,int(math.sqrt(nvariables))):
					pep_v_new.append(matriz2[i][j])
					
			# 7) Calculate determinants of nugget effect and check if LMC is found			

			determ2 = np.linalg.det(matriz2_new)
		
			if (determ2 > 0):
				check_LMC2 = True
			
			count = count + 1
			if (count == 200):
				print("LMC not found")
				check_LMC = True
				check_LMC2 = True
		#...................................................................................

		'''
		ADOPT MODEL OF PRINCIPAL CONTINUITY FOR SECONDARY AND VERTICAL CONTINUITY 
		_________________________________________________________________________

		This part of algorithm normalizes the secondary and vertical directions with the 
		same parameters of principal continuity fitting. The following steps are:

		1) Calculate average range for all variables 
		2) Adopt principal continuity model for secondary and vertical directions 
		3) Adopt the same contribution and nugget effect of principal direction
		4) Print on command pannel the variogram parameters 	
		5) Plot variograms in pannel	

		'''



		# 1) Calculate average range for all variables 
	
		
		rangemedio = [0 for i in range(0,len(rangeV[0]))]		
		for m in rangeV:				
			for i in range(0,len(m)):
				rangemedio[i] = rangemedio[i] + m[i]/len(rangeV)
			
		
			
		# 2) Adopt principal continuity model for secondary and vertical directions 			
			
		modelo_adotado = []			
		modelo_adotado = testerV[0] 
			
		for i in range(0, len(testerV)):
			testerV[i] = modelo_adotado

		# 3) Adopt the same contribution and nugget effect of principal direction

		if (direcao > 1):
			pepV = pepm
			testerV = modeling
			sillV = sillm
			
		# 4) Print on command pannel the variogram parameters 

		if (direcao ==1):	
			print("Variogram parameters")
			print(".............................................................................")
			print("model = 0 -> spherical, model = 1 -> exponencial, model =2 -> gaussian")
			print(".............................................................................")					
			print("sills:" + str (matriz_contr)) 
			print("models:" + str(testerV[0]))
			print("nugget:" + str(matriz2_new))
			print("range: da direcao "+str(direcao) +"  "+ str(rangemedio))
		else:
			print("range: da direcao "+str(direcao) +"  "+ str(rangemedio))
			print(".............................................................................")
						
	
		# 5) Plot variograms in pannel
	
		for n in range(0,nvariables):
			r_sill = sillV[n]
			r_range = rangemedio[0]
			r_pep = pep_v_new[n]
			r_testr = testerV[n]
			variogramas = variogramasV[n]
			lags = lagsV[n]
			variograma = variogramas[I_dir_principal]
			lag = lags[I_dir_principal]

			funcao = []
			for i in np.linspace(0, max(lag), 400):
		    		valor_variograma = 0
		    		for j in range(0, numero_estruturas):
		        		if (r_testr[j] == 0):
		            			if (i <= r_range):
		                			valor_variograma = valor_variograma + r_sill[j]*(1.5*i/r_range-0.5*math.pow((i/r_range),3))
		            			else:
		                			valor_variograma = valor_variograma + r_sill[j]
		        		elif (r_testr[j] == 1):
		            			valor_variograma = valor_variograma + r_sill[j]*(1-math.exp(-3*i/r_range))
		        		elif (r_testr[j] == 2):
		            			valor_variograma = valor_variograma + r_sill[j]*(1-math.exp(-3*pow(i/r_range,2)))
		    		valor_variograma = valor_variograma + r_pep
		    		funcao.append(valor_variograma)
			limite = np.linspace(0, max(lag),400)
			
			plt.subplot(math.sqrt(nvariables),math.sqrt(nvariables),n+1)			
			plt.plot(lag,variograma,'g^--',limite, funcao)
			
			if (direcao == 1):
				plt.title( "\n \n Principal"+str(indice_head[n]) +" " +str(indice_tail[n]))	
			if (direcao == 2):
				plt.title( "\n \n Secondary"+str(indice_head[n]) +" " +str(indice_tail[n]))
			if (direcao == 3):
				plt.title( "\n \n Vertical" +str(indice_head[n]) +" " +str(indice_tail[n]))
		plt.show()	
		return testerV, pepV, sillV_new	
		# ........................................................................................

	# Call variogram optimization for each direction
		
	if (set_maximumrange == 1 and set_minrange == 0 and set_verticalrange == 0):		
		pepV =[]		
		saida, pepV, sillV = otimizacao(I_dir_principal,variogramasV, lagsV, paresV, nvariables, numero_estruturas, restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range_max,max_range_max, -1, 1, -1, [] )
	if (set_maximumrange == 1 and set_minrange == 1 and set_verticalrange == 0):
		pepV =[]
		saida, pepV, sillV = otimizacao(I_dir_principal,variogramasV, lagsV, paresV, nvariables, numero_estruturas, restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range_max,max_range_max, -1, 1, -1, [] )
		saida2, pepd, sillV = otimizacao(I_dir_secundario,variogramasV, lagsV, paresV, nvariables, numero_estruturas, restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range_min,max_range_min, saida, 2, pepV, sillV )
	if (set_verticalrange == 1 and set_maximumrange ==1 and set_minrange == 1 ):  	
		pepV =[]
		saida, pepV, sillV = otimizacao(I_dir_principal,variogramasV, lagsV, paresV, nvariables, numero_estruturas, restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range_max,max_range_max, -1, 1, -1, [] )
		saida2, pepd, sillV = otimizacao(I_dir_secundario,variogramasV, lagsV, paresV, nvariables, numero_estruturas, restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range_min,max_range_min, saida, 2, pepV, sillV )	
		saida3, pepe, sillV = otimizacao(I_dir_vertical,variogramasV, lagsV, paresV, nvariables, numero_estruturas, restriction,min_contribution,max_contribution,min_nugget,max_nugget,min_range_vert,max_range_vert, saida, 3, pepV, sillV )     		

    def finalize(self):
        print "Finalize optimizing_variogram"
        return True
    def name(self):
        return "automatic_fitting2"
##############################################
def get_plugins():
    return["automatic_fitting2"]
