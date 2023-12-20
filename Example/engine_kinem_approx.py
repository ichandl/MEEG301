# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 19:03:12 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# geometric parameters
w = 1500/60*2*np.pi                  # rotation rate [rad/s]
t = np.linspace(0,4*np.pi/w,101)     # time array [s]
R2 = 0.05                            # crank length [m]
R3 = 0.20                            # connecting rod length [m]

# complex positions
tht3 = np.arcsin(-R2/R3*np.sin(w*t))     # connecting rod angle [rad]
z_A = R2*np.exp(1j*w*t)                  # position of point A [m]
z_B = z_A + R3*np.exp(1j*tht3)           # position of point B [m]

# plot trajectory
plt.figure()
plt.plot(np.real(z_A),np.imag(z_A),label='crank')
plt.plot(np.real(z_B),np.imag(z_B),label='piston')
plt.axis('equal')
plt.legend()

# piston position spline curve fit
x_B_spline = interpolate.splrep(t,np.real(z_B))   # spline of piston x-position array

# approximate piston position and its spline
x_B_approx = R3 - R2**2/(4*R3) + R2*(np.cos(w*t) + R2/(4*R3)*np.cos(2*w*t))
x_B_approx_spline = interpolate.splrep(t,x_B_approx)   # spline of piston x-position array

# plot comparison between x_B and its approximation
plt.figure()
plt.plot(t,np.real(z_B)*100,label='exact')
plt.plot(t,x_B_approx*100,label='approx')
plt.ylabel('Point B x position [cm]')
plt.xlabel('Time [s]')
plt.legend()

plt.figure()
plt.plot(t,(np.real(z_B)-x_B_approx)*1000)
plt.ylabel('Approximation error [mm]')
plt.xlabel('Time [s]')

plt.show()