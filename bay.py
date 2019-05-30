#
#	box model of the carbon cycle in a system composed by a bay, surface and deep global oceans
#
#       boxes are as follows, 1=BAY, 2=SGO, 3=DGO
#
#       agostino merico
#       started: 02 march 2013
#
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# riverine input

RC=0.001 

# fluxes between reservoirs

SB=0.05  
BS=0.02
SD=0.40
DS=20.0

# set of ODEs

def func(y, t):

    tc=3.15576e-4         # converts from s-1 to yr-1
    V=[0.0001, 1.5, 12.0] # volumes of the 3 reservoirs in 1e17 m3

    dydt=np.zeros(3)
    
    dydt[0]=((RC+SB*y[1]-BS*y[0])/V[0])*tc              # BAY
    dydt[1]=((BS*y[0]-SB*y[1]+DS*y[2]-SD*y[1])/V[1])*tc # SOG
    dydt[2]=((SD*y[1]-DS*y[2])/V[2])*tc                 # DOG

    return dydt

yi=[0.01, 10., 20.]            # initial concentrations in the three reservoirs (in mmol m-3)
t=np.linspace(0., 1500., 100.) # time grid

# solve ODEs

sol=odeint(func,yi,t)

BAY=sol[:,0]
SOG=sol[:,1]
DOG=sol[:,2]

# plot results

plt.figure()

plt.plot(t,BAY, label='bay')
plt.plot(t,SOG, label='surface global ocean')
plt.plot(t,DOG, label='deep global ocean')

plt.xlabel('time (years)')
plt.ylabel('carbon (mmol m$^{-3}$)')
plt.legend(loc=0)


plt.savefig('c3box.pdf')

plt.show()
