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
	return ((j*m*r / (4*math.pi)) *(-1) * (np.sin(p)*(y-w)) / ( ( ((x)-r * np.cos(p))**2 +((-z)-r*np.sin(p))**2 + (y-w)**2 )**(1.5) ))

def fny(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi))  * (r-(((-z)*np.sin(p))+((x)*np.cos(p)))) / ( ( ((x)-r * np.cos(p))**2 +((-z)-r*np.sin(p))**2 + (y-w)**2 )**(1.5) ))


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
#Konstanten und Integrations und Plottingparameter
m = 1.256637061e-6
#m=1
length=1

aw=-d/2 #Integrationsgrenzen fuer Volumenintegral
bw= d/2
ap=0.0 #
bp=2*math.pi #
ar=ri #
br=ra #

def f1(abc):
	return ap

def f2(abc):
	return bp

def f3(abc, cde):
	return ar

def f4(abc, cde):
	return br

#Multiprocessing init Lists
list1=(fnx, aw, bw, f1, f2, f3, f4, (x, y, z, j, m)) #1e-05, 1e-03)
list2=(fny, aw, bw, f1, f2, f3, f4, (x, y, z, j, m)) #1e-05, 1e-03)
list3=(fnz, aw, bw, f1, f2, f3, f4, (x, y, z, j, m)) #1e-05, 1e-03)

t1=time.time()
for s in range(5, 6):
	x=s/50
	print("----------- S = ", s)
	for t in range(-11, -5):
		y=t/50
		print("- T = ", t)
		pool=mp.Pool(3)
		results=pool.starmap(tplquad, [(fnx, aw, bw, f1, f2, f3, f4, (x, y, z, j, m), 4.49e-06, 1.49e-04), (fny, aw, bw, f1, f2, f3, f4, (x, y, z, j, m), 4.49e-06, 1.49e-04), (fnz, aw, bw, f1, f2, f3, f4, (x, y, z, j, m), 4.49e-06, 1.49e-04)])
		pool.close()
		pool.join()
		Ix=results[0]
		Iy=results[1]
		Iz=results[2]
		
		#Ix=tplquad(fnx, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m))# epsabs=1.0e-06, epsrel = 1.0e-05)
		#Iy=tplquad(fny, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m))# epsabs=1.0e-06, epsrel = 1.0e-05)
		#Iz=tplquad(fnz, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, (x, y, z, j, m))# epsabs=1.0e-06, epsrel = 1.0e-05)
		abs = math.sqrt(Ix[0]**2 + Iy[0]**2 + Iz[0]**2)   #eventuell noch Begrenzung nach oben/unten
		#print(x, y, z, Ix[0], Iy[0], Iz[0], abs)
		tf = open( 'datafile2.txt', 'a')
		tf.write(repr(x)+" " +repr(y) +" "+repr(z) + " "+ repr(Ix[0]*length/abs) + " " + repr(Iy[0]*length/abs) + " " + repr(Iz[0]*length/abs)+ " " + repr(140*abs) +'\n') #3D
		tf.close()
t2=time.time()
print(t2-t1)
