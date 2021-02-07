# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 09:46:42 2021

@author: RBene
"""
#Program for SE 262 Homework 2: Stepped Lap Joints
#Derived using the Volkersen Model
#Calculates the Internal force Profile (Nox) over Length of Adhesive
import math
import numpy as np
import matplotlib.pyplot as plt
n=4 #number of steps
m=n-1 #number of steps, make sure to -1 for python indexing
Ga=180*1000 #shear modulus adhesive
ta=.006 #thickness adhesive
Eo=11.5*1000000 #outer adherend Y modulus
Ei=10*1000000 #inner adherend Y modulus
to=[.008,.016,.024,.032] #outer adherend thickness
ti=[.036,.036,.036,.036] #inner adherend thickness
c=[.25,.25,.25,.9] #overlap per each lap
Nx=500 #tensile applied load
deltaT=0 #change in temp
i=0
lamda=[]
Cok=[]
alpha=[]
beta=[]
hi=[]
No=[]
Aok=[]
Bok=[]
Nox=[]
j=0
p=0
#creates function that iterates through and gets your Nox
while p<n:
    while j<m:
        while i<(len(ti)):
            lamda.append(math.sqrt((Ga/ta)*((1/(Eo*to[i]))+(1/(Ei*ti[i]))))) #calculate lambda
            Cok.append(Ga*Nx/(Ei*ta*ti[i]))
            alpha.append(lamda[i]*math.tanh(lamda[i]*c[i])/2)
            beta.append(lamda[i]/(2*math.tanh(lamda[i]*c[i])))
            hi.append(2*alpha[i]*Cok[i]/lamda[i]**2)
            i+=1
            #print('If this values is less than 3, simplify matrix=', lamda[i]*c)
        #Can only perform this step below if simplifed matrix above
        Nox.append((hi[j]+hi[(j+1)])/(lamda[j]+lamda[(j+1)])) #adjust per given m value
        j+=1
        #Boundary conditions for Nox
    if j==m:
        Nox.insert(0,0) #side without any force
        Nox.insert(len(Nox),Nx) #side with applied force
    j+=1
    x=np.arange(-c[p],(c[p]),0.01) #x values
    for t in x:    
        Aok=(Nox[p]+Nox[(p+1)]-(2*Cok[p]/lamda[p]**2))/(2*math.cosh(lamda[p]*c[p]))
        Bok=(Nox[(p+1)]-Nox[p])/(2*math.sinh(lamda[p]*c[p]))
        No.append(Aok*math.cosh(lamda[p]*t)+Bok*math.sinh(lamda[p]*t)+Cok[p]/lamda[p]**2)
    p+=1

#Set up iteration for graphing
x1=np.arange(0,(sum(c)*2),0.01)  #will need to change for larger steps

#plotting
plt.plot(x1,No,'r',label='Stepped Lap Joint')
plt.xlabel('Location on Adhesive (in.)')
plt.ylabel('No (lbf/in)')
plt.title("Internal Force Profile over Length of Adhesive")
plt.legend(loc="best")
plt.grid(True)
