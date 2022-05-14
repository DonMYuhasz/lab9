# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 22:08:59 2022

@author: DOnald Yuhasz
Modified from code on http://apmonitor.com/che263/index.php/Main/PythonDynamicSim
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# start values
Initial_Susceptible = 100
Initial_Infected = 10

#short function to prevent negative numbers of people
def nozero(i,y):
    if y+i < 0:
        return -y
    return i

def sidarthe(m,t):
   
   # assign array to var for readability
   S = m[0]
   I = m[1]
   D = m[2]
   A = m[3]
   R = m[4]
   T = m[5]
   H = m[5]
   E1 = m[5]
   #transfer coefficients
   
   a = .01#alpha
   B = .01 #beta
   y1 = .01 #rho?
   d = .01 #delta (lowercase)
   e = .01 #epsilon
   s = .01 #I genuinly don't know
   v = .01 #gamma
   n = .01
   p = .01
   theta = .01
   u = .01 #mu
   k = .01
   V = .01
   E = 0.01
   t = 0.1
   o = 0.01
   
   
   # differential equations
   
   susceptible = -S*((a*I)+(B*D)+(y1*A)+(d*R))
   infected = S*((a*I)+(B*D)+(y1*A)+(d*R)) - (I*(e +s+ v))
   diagnosed = e* I - (D*(n+p)) 
   ailing = s*I - (A*(theta + u + k))
   recognized = (n*D) + (theta * A)  - (R *(V+E))
   thretened = (u*A)+(v*R) -(T*(o + t))
   healed = (v * I) + (p*D) + (k*A) + (E*R)+ (o*T)
   extinct = t*T
   
   #prevents negative people
   susceptible = nozero(susceptible,S)
   infected = nozero(infected,I)
   diagnosed = nozero(diagnosed,D)
   ailing  = nozero(ailing,A)
   recognized = nozero(recognized,R)
   thretened = nozero(thretened,T)
   healed = nozero(healed,H)
   extinct = nozero(extinct,E1)

   # no overflow conditions
   #should not be possible to underflow 
   
   dhdt = [susceptible,infected,diagnosed,ailing,recognized,thretened,healed,extinct]
   return dhdt

# integrate the equations
t = np.linspace(0,100) # times to report solution
P0 = [Initial_Susceptible,Initial_Infected,0,0,0,0,0,0]            # initial conditions for no. people
y = odeint(sidarthe,P0,t) # integrate

# plot results
plt.figure(1)
plt.plot(t,y[:,0],'b-')
plt.plot(t,y[:,1],'r--')
plt.plot(t,y[:,2],'g')
plt.plot(t,y[:,3],'p')
plt.plot(t,y[:,4],'g')
plt.plot(t,y[:,5],'g')
plt.legend(['susceptible','infected','diagnosed','ailing','recognized'])
plt.figure(2)
plt.plot(t,y[:,0],'b')
plt.plot(t,y[:,6],'g')
plt.plot(t,y[:,7],'r')
plt.xlabel('Time (hrs)')
plt.ylabel('occupants')
plt.legend(['susceptible','healed','extinct'])
plt.show()