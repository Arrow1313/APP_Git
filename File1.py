__author__ = 'Gruppe_B'
#Numerische Integration

from scipy.integrate import quad
import numpy as np
import math


def fnx(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (np.cos(p)*z) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

def fny(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (np.sin(p)*z) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

def fnz(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (r-((y*np.sin(p))+(x*np.cos(p)))) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

x=0.0
y=0.0
z=0.0
r=0.09
i=2.0
m=1.256637061e-6
a=0.0
b= 2*math.pi

for l in range(-3, 4):__author__ = 'Gruppe_B'
#Numerische Integration

from scipy.integrate import quad
import numpy as np
import math


def fnx(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (np.cos(p)*z) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

def fny(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (np.sin(p)*z) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

def fnz(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (r-((y*np.sin(p))+(x*np.cos(p)))) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

x=0.0
y=0.0
z=0.0
r=0.09
i=2.0
m=1.256637061e-6
a=0.0
b= 2*math.pi__author__ = 'Gruppe_B'
#Numerische Integration

from scipy.integrate import quad
import numpy as np
import math


def fnx(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (np.cos(p)*z) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

def fny(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (np.sin(p)*z) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

def fnz(p,x,y,z,r,i,m):
    return ((i*m*r / (4*math.pi)) * (r-((y*np.sin(p))+(x*np.cos(p)))) / ( ( (x-r * np.cos(p))**2 +(y-r*np.sin(p))**2 + z**2 )**(1.5) ))

x=0.0
y=0.0
z=0.0
r=0.09
i=2.0
m=1.256637061e-6
a=0.0
b= 2*math.pi

for l in range(-3, 4):
    z=l/10
    for k in range(-3,4):
        x=k/10
        Ix=quad(fnx, a, b, args=(x, y, z, r, i, m))
        Iy=quad(fny, a, b, args=(x, y, z, r, i, m))
        Iz=quad(fny, a, b, args=(x, y, z, r, i, m))
        print("[", Ix[0], ",", Iy[0], ",", Iz[0], "]")

for l in range(-3, 4):
    z=l/10
    for k in range(-3,4):
        x=k/10
        Ix=quad(fnx, a, b, args=(x, y, z, r, i, m))
        Iy=quad(fny, a, b, args=(x, y, z, r, i, m))
        Iz=quad(fny, a, b, args=(x, y, z, r, i, m))
        print("[", Ix[0], ",", Iy[0], ",", Iz[0], "]")
    z=l/10
    for k in range(-3,4):
        x=k/10
        Ix=quad(fnx, a, b, args=(x, y, z, r, i, m))
        Iy=quad(fny, a, b, args=(x, y, z, r, i, m))
        Iz=quad(fnz, a, b, args=(x, y, z, r, i, m))
        print("[", Ix[0], ",", Iy[0], ",", Iz[0], "]")