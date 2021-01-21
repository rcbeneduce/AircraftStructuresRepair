#Author: Ryan Beneduce
#For SE 264 UCSD
import math
import numpy as np
import matplotlib.pyplot as plt

#Aircraft structural repair
#Volkersen Model for Lapped Joints
Ga=212*1000 #shear modulus for adhesive
Eo=12.5*1000000 #Young's Mod for outer adherend
Ei=17*1000000 #Young's Mod for inner adherend
ti=.080 #thickness inner adherand
ta=.010 #thickness adhesive
to=.075 #thermal application
alphao=1.1*10**-6 #thermal expansion for outer adherand
alphai=5*10**-6 #thermal expansion for inner adherand
deltaT=300-70 #change in temp
lamda=math.sqrt((Ga/ta)*((1/(Eo*to))+(1/(Ei*ti)))) #calculate lambda
shear=15000 #ultimate shear stress
#calculate max Nx
c=0.5 #range of model
x=np.arange(-c,c+.01,0.01) #x values
shearx=[]
for t in x:
    Nx=500
    Co=(Ga/ta)*((Nx)/(Ei*ti)-deltaT*(alphao-alphai))
    a=(math.sinh(lamda*t))/(math.cosh(lamda*c))
    b=(Nx/2)*(math.cosh(lamda*t))/(math.sinh(lamda*c))
    d=Nx/2-Co/lamda**2
    sheary=lamda*((d*a)+b) #shearx based on x
    shearx.append(sheary)
#type C
to=0.015
lamda=math.sqrt((Ga/ta)*((1/(Eo*to))+(1/(Ei*ti))))
shearx2=[]
for w in x:
    p=math.sinh(lamda*w)/(2*math.cosh(lamda*c))
    q=(Ga*math.sinh(lamda*w))/(Ei*ta*ti*lamda**2*math.cosh(lamda*c))
    r=(math.cosh(lamda*w))/(2*math.sinh(lamda*c))
    yes=700
    sheary=lamda*(yes*p-yes*q+yes*r) #shearx based on x
    shearx2.append(sheary)
    
plt.plot(x,shearx,'r')
plt.plot(x,shearx2,'b',label="Type C")
plt.xlabel('X(in.)')
plt.ylabel('Adhesive Shear Stress (psi)')
plt.title("Adhesive Shear Stress Profile with ")
plt.grid(True)
