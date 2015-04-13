__author__ = 'Gruppe_B'
#Numerische Dreifachintegration

from scipy.integrate import tplquad
import numpy as np
import math
from functools import partial
import time
import multiprocessing as mp


#Integrale
def fnx(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi)) * (np.cos(p)*(y-w)) / ( ( ((x)-r * np.cos(p))**2 +((-z)-r*np.sin(p))**2 + (y-w)**2 )**(1.5) ))

def fnz(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi)) * (np.sin(p)*(y-w)) / ( ( ((x)-r * np.cos(p))**2 +((-z)-r*np.sin(p))**2 + (y-w)**2 )**(1.5) ))

def fny(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi)) * (-1) * (r-(((-z)*np.sin(p))+((x)*np.cos(p)))) / ( ( ((x)-r * np.cos(p))**2 +((-z)-r*np.sin(p))**2 + (y-w)**2 )**(1.5) ))


#Koordinaten
x=0.0
y=0.0
z=0.0

#Spuleneigenschaften
ri=0.1
dr=0.01
ra=ri+dr #
d=0.32
i=1
n=200
j=i*n/(d*dr) #
#Konstanten und Integrations- und Plottingparameter
m = 1.256637061e-6
#m=1
length=1

aw=-d/2 #Integrationsgrenzen für Volumenintegral
bw= d/2
ap=0.0 #
bp=2*math.pi #
ar=ri #
br=ra #

#Multiprocessing init Lists
#lista=[fnx, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m)]
#listb=[fny, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m)]
#listc=[fnz, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m)]

t1=time.time()
for s in range(10, 12):
	x=s/50
	print("----------- S = ", s)
	for t in range(-11, 12):
		y=t/50
		print("- T = ", t)
		#pool=mp.Pool(1)
		#results=pool.map(tplquad, [fnx, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m)])
		#pool.close()
		#pool.join()
		#Ix=results[0]
		#Iy=results[1]
		#Iz=results[2]
		Ix=tplquad(fnx, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m), epsabs=1.0e-05, epsrel = 1.0e-04)
		Iy=tplquad(fny, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m), epsabs=1.0e-05, epsrel = 1.0e-04)
		Iz=tplquad(fnz, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m), epsabs=1.0e-05, epsrel = 1.0e-04)
		abs = math.sqrt(Ix[0]**2 + Iy[0]**2 + Iz[0]**2)   #eventuell noch Begrenzung nach oben/unten
		#print(x, y, z, Ix[0], Iy[0], Iz[0], abs)
		print(Ix, Iy)
		#tf = open( 'datafile2.txt', 'a')
		#tf.write(repr(x)+" " +repr(y) +" "+repr(z) + " "+ repr(Ix[0]*length/abs) + " " + repr(Iy[0]*length/abs) + " " + repr(Iz[0]*length/abs)+ " " + repr(140*abs) +'\n') #3D
		#tf.close()
t2=time.time()
print(t2-t1)
