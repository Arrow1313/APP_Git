__author__ = 'Gruppe_B'
#Numerische Dreifachintegration

from scipy.integrate import tplquad
import numpy as np
import math

#Integrale
def fnx(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi)) * (np.cos(p)*(-y-w)) / ( ( ((x)-r * np.cos(p))**2 +((z)-r*np.sin(p))**2 + (-y-w)**2 )**(1.5) ))

def fnz(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi)) *(-1) * (np.sin(p)*(-y-w)) / ( ( ((x)-r * np.cos(p))**2 +((z)-r*np.sin(p))**2 + (-y-w)**2 )**(1.5) ))

def fny(r,p,w,x,y,z,j,m):
	return ((j*m*r / (4*math.pi)) * (r-(((z)*np.sin(p))+((x)*np.cos(p)))) / ( ( ((x)-r * np.cos(p))**2 +((z)-r*np.sin(p))**2 + (-y-w)**2 )**(1.5) ))


#Koordinaten
x=0.0
y=0.0
z=0.0

#Spuleneigenschaften
ri=0.1
dr=0.01
ra=ri+dr #
d=0.2
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

for s in range(-10, 11):
	x=s/50
	print(s)
	for t in range(-15, 16):
		y=t/50
		Ix=tplquad(fnx, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, args=(x, y, z, j, m))
		Iy=tplquad(fny, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, args=(x, y, z, j, m))
		Iz=tplquad(fnz, aw, bw, lambda w: ap, lambda w: bp, lambda w,p: ar, lambda w,p: br, args=(x, y, z, j, m))
		abs = math.sqrt(Ix[0]**2 + Iy[0]**2 + Iz[0]**2)   #eventuell noch Begrenzung nach oben/unten
		#print(x, y, z, Ix[0], Iy[0], Iz[0], abs)
		tf = open( 'datafile2.txt', 'a')
		tf.write(repr(x)+" " +repr(y) +" "+repr(z) + " "+ repr(Ix[0]*length/abs) + " " + repr(Iy[0]*length/abs) + " " + repr(Iz[0]*length/abs)+ " " + repr(140*abs) +'\n') #3D
		tf.close()