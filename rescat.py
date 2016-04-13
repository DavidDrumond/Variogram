#!/bin/python
import sgems
import numpy as np
import math
import random




"""
calculating distace
"""

class rescat:
    def _init_(self):
        pass 
    def initialize(self,params):
        self.params = params
        return True
    def execute(self):
	print self.params
	# ESTABELECA OS PARAMETROS INICIAIS
	propertie = self.params['prop']['property']
	grid = self.params['prop']['grid']
	
	measured = float(self.params['Measured']['value'])
	Indicated = float(self.params['Indicated']['value'])
	Inferred = float(self.params['Inferred']['value'])
	dx = float(self.params['dx']['value'])
	dy = float(self.params['dy']['value'])
	dz = float(self.params['dz']['value'])
	
	grid = self.params['prop']['grid'] 
	x = sgems.get_property(grid, "_X_")
	y = sgems.get_property(grid, "_Y_")
	z = sgems.get_property(grid, "_Z_")
	distance = []
	
	

	numx = int(round((max(x) - min(x))/dx,0)) + 1
	numy = int(round((max(y) - min(y))/dy,0)) + 1
	numz = int(round((max(z) - min(z))/dz,0)) + 1

		

	cmd="NewCartesianGrid rescat::"+str(numx)+"::"+str(numy)+"::"+str(numz)+"::"+str(dx)+"::"+str(dy)+"::"+str(dz)+"::"+str(min(x))+"::"+str(min(y))+"::"+str(min(z))+"::0"
        print cmd
	sgems.execute(cmd)

	xh = sgems.get_property("rescat","_X_")
	yh = sgems.get_property("rescat","_Y_")
	zh = sgems.get_property("rescat","_Z_")	

	

	
	v = []
	for i, j, k in zip(xh,yh,zh):
		valores = 0
		menor_valor = 1000
		for xind, yind, zind in zip(x,y,z):
			distance = (math.sqrt((xind-i)**2+(yind-j)**2+(zind-k)**2))
			if (distance < measured):
				valores = 0
			elif( distance < Indicated and distance > measured):
				valores = 1
			elif(distance < Inferred and distance > Indicated):
				valores = 2
			else:
				valores = 3
			if (valores < menor_valor):
				menor_valor = valores		
		v.append(menor_valor)
	print (v)
	sgems.set_property("rescat","R",v)
	
	
    	return True 

    def finalize(self):
        print "Finalize program"
        return True
    def name(self):
        return "rescat"
##############################################
def get_plugins():
    return["rescat"]
