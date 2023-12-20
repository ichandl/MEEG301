# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 11:15:43 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate


w = np.pi                          # rotation rate [rad/s]
t = np.linspace(0,2*np.pi/w,101)   # time array [s]
R = 0.1                            # radius [m]
z = R*np.exp(1j*w*t)               # complex position [m]

# plot trajectory
plt.figure()
plt.plot(np.real(z),np.imag(z))
plt.axis('equal')

# analytical velocity
dzdt = R*1j*w*np.exp(1j*w*t)       # complex velocity [m/s]



# plot x-y position and velocity vs time
_, ax = plt.subplots(2,1)

ax[0].plot(t,np.real(z),label='x-position')
ax[0].plot(t,np.real(dzdt),label='x-velocity')
ax[0].set_xlabel('Time')
ax[0].legend()

ax[1].plot(t,np.imag(z),label='y-position')
ax[1].plot(t,np.imag(dzdt),label='y-velocity')
ax[1].set_xlabel('Time')
ax[1].legend()



# quiver plot
t_q = np.linspace(0,2*np.pi/w,9)   # coarser time array [s]
z_q = R*np.exp(1j*w*t_q)           # complex position [m]
dzdt_q = R*1j*w*np.exp(1j*w*t_q)   # complex velocity [m/s]
plt.figure()
plt.plot(np.real(z),np.imag(z))
plt.quiver(np.real(z_q),np.imag(z_q),np.real(dzdt_q),np.imag(dzdt_q))
plt.axis('equal')

