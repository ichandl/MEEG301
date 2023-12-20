# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 19:03:12 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# crank position
w = np.pi                            # rotation rate [rad/s]
t = np.linspace(0,4*np.pi/w,101)     # time array [s]
R2 = 5                               # radius [cm]
z_A = R2*np.exp(1j*w*t)              # complex position [cm]

# slider position
R3 = 15                                     # coupler length [cm]
D = 8                                       # slider vertical offset [cm]
tht3 = np.arcsin((D-R2*np.sin(w*t))/R3)     # coupler angle [rad]
z_B = z_A + R3*np.exp(1j*tht3)              # complex position [cm]

# plot trajectory
plt.figure()
plt.plot(np.real(z_A),np.imag(z_A),label='crank')
plt.plot(np.real(z_B),np.imag(z_B),label='slider')
plt.axis('equal')
plt.legend()

# spline curve fits
x_A_spline = interpolate.splrep(t,np.real(z_A))   # spline of crank x-position array
y_A_spline = interpolate.splrep(t,np.imag(z_A))   # spline of crank y-position array
x_B_spline = interpolate.splrep(t,np.real(z_B))   # spline of slider x-position array
y_B_spline = interpolate.splrep(t,np.imag(z_B))   # spline of slider y-position array

# numerically differentiate spline fits
_, ax = plt.subplots(2,1)

ax[0].plot(t,interpolate.splev(t,x_A_spline,der=1),label='crank x-velocity')
ax[0].plot(t,interpolate.splev(t,x_B_spline,der=1),label='slider x-velocity')
ax[0].set_xlabel('Time')
ax[0].legend()

ax[1].plot(t,interpolate.splev(t,y_A_spline,der=1),label='crank y-velocity')
ax[1].plot(t,interpolate.splev(t,y_B_spline,der=1),label='slider y-velocity')
ax[1].set_xlabel('Time')
ax[1].legend()