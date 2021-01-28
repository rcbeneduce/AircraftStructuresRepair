-Part A (Axial force in outer adherand)
import math
import numpy as np
import matplotlib.pyplot as plt

Ga=112*1000
Eo=3.7*1000000
Ei=3.7*1000000
ti=.062
ta=.004
to=.015
Nx=1660
c=.25
t=np.arange(-c,c+.01,.01)
lamda=math.sqrt((Ga/ta)*((1/(Eo*to))+(1/(Ei*ti))))
Co=(Ga*Nx)/(Ei*ta*ti)
Ao=(Nx-(Co/lamda**2))*(1/math.cosh(lamda*c))
N=[]
for x in t:
    No=Ao*math.cosh(lamda*x)+Co/lamda**2
    N.append(No)
plt.plot(t,N,'r')
plt.grid(True)


-Part B (Shear stress)

import math
import numpy as np
import matplotlib.pyplot as plt

Ga=112*1000
Eo=3.7*1000000
Ei=3.7*1000000
ti=.062
ta=.004
to=.015
Nx=1660
c=.25
t=np.arange(-c,c+.01,.01)
lamda=math.sqrt((Ga/ta)*((1/(Eo*to))+(1/(Ei*ti))))
Co=(Ga*Nx)/(Ei*ta*ti)
Ao=(Nx-(Co/lamda**2))*(1/math.cosh(lamda*c))
N=[]
for x in t:
    shearstress=lamda*Ao*math.sinh(lamda*x)
    N.append(shearstress)
plt.plot(t,N,'r')
plt.grid(True)
