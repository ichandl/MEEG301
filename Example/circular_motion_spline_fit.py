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
a = 15*np.pi/180                   # amplitude of rocker motion [rad]
b = 70*np.pi/180                   # offset of rocker motion [rad]
z = R*np.exp(1j*(a*np.sin(w*t)+b)) # complex position [m]

# plot trajectory
plt.figure()
plt.plot(np.real(z),np.imag(z))
plt.axis('equal')

# analytical velocity
dzdt = R*1j*a*w*np.cos(w*t)*np.exp(1j*(a*np.sin(w*t)+b))       # complex velocity [m/s]

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


# spline curve fits
x_spline = interpolate.splrep(t,np.real(z))   # spline of x-position array
y_spline = interpolate.splrep(t,np.imag(z))   # spline of y-position array

# evaluate spline fits and compare to original position data
_, ax = plt.subplots(2,1)

ax[0].plot(t,np.real(z),label='x-position')
ax[0].plot(t,interpolate.splev(t,x_spline),'r--',label='x-spline')
ax[0].set_xlabel('Time')
ax[0].legend()

ax[1].plot(t,np.imag(z),label='y-position')
ax[1].plot(t,interpolate.splev(t,y_spline),'r--',label='y-spline')
ax[1].set_xlabel('Time')
ax[1].legend()

# numerically differentiate spline fits and compare to original velocity data
_, ax = plt.subplots(2,1)

ax[0].plot(t,np.real(dzdt),label='x-velocity')
ax[0].plot(t,interpolate.splev(t,x_spline,der=1),'r--',label='dx-spline')
ax[0].set_xlabel('Time')
ax[0].legend()

ax[1].plot(t,np.imag(dzdt),label='y-velocity')
ax[1].plot(t,interpolate.splev(t,y_spline,der=1),'r--',label='dy-spline')
ax[1].set_xlabel('Time')
ax[1].legend()
