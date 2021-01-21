# -*- coding: utf-8 -*-
"""
Created on Sat Jan 16 17:19:28 2021

@author: Ryan
"""
import math
import numpy as np
import matplotlib.pyplot as plt
Ga=212*1000
Eo=12.5*1000000
Ei=17*1000000
ti=.080
ta=.010
to=.075 #thermal application
alphao=1.1*10**-6
alphai=5*10**-6
deltaT=300-70
lamda=math.sqrt((Ga/ta)*((1/(Eo*to))+(1/(Ei*ti)))) #calculate lambda
shear=5500 #ultimate shear
#calculate max Nx
c=0.5
x=np.arange(-c,c+.01,0.01) #x values
Nx=[]
for t in x:
    lamda=math.sqrt((Ga/ta)*((1/(Eo*to))+1/(Ei*ti)))
    a=(math.sinh(lamda*t))/(2*math.cosh(lamda*c))
    b=(Ga*math.sinh(lamda*t))/(Ei*ta*ti*lamda**2*math.cosh(lamda*c))
    d=(math.cosh(lamda*t))/(2*math.sinh(lamda*c))
    e=Ga*deltaT*(alphao-alphai)*math.sinh(lamda*t)
    f=lamda**2*math.cosh(lamda*c)*ta
    Nx1=-1*((shear/lamda)-(e/f))/(a-b+d)
    Nx.append(Nx1)

plt.plot(x,Nx,'r')
plt.xlabel('Location on Adhesive (in.)')
plt.ylabel('Joint Failure Load (lbs/in)')
plt.title("Variable Failure Load vs Location on Adhesive")
#plt.legend(loc="best")
plt.grid(True)
