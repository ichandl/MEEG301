# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 22:18:20 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy import interpolate

# givens
R1 = 10                              # ground link length [cm]
R2 = 4                               # crank length [cm]
R3 = 12                              # coupler length [cm]
R4 = 8                               # rocker length [cm]
w_f = 2*np.pi                        # final (steady state) angular rate [rad/s]
tau = 1.5                            # crank time constant [s]
link_config = 'open'                 # 'open' or 'crossed'
z1 = R1                              # ground link

# set up arrays for simulation
t = np.linspace(0,10,5001)            # time array [s]
w = w_f*(1-np.exp(-t/tau))           # crank angular rate array [rad/s]
tht2 = np.zeros_like(t)              # crank angle array [rad]
tht3 = np.zeros_like(t)              # coupler angle array [rad]
tht4 = np.zeros_like(t)              # rocker angle array [rad]

plt.figure()
plt.plot(t,w)
plt.ylabel('Crank rotation rate [rad/s]')
plt.xlabel('Time [s]')

# choose which solution of tht4 to look for
if link_config == 'open':
    bracket = [0,np.pi]
    x0 = np.pi/2
else:
    bracket = [-np.pi,0]
    x0 = -np.pi/2


# function whose root we want to find within the bracket
def calc_tht4(tht4_guess):
    z4 = R4*np.exp(1j*tht4_guess)
    z5 = z1 + z4
    z3 = z5 - z2
    r = np.abs(z3)
    return r - R3


for i in range(t.size):
    if i == 0:
        tht2[i] = 0
    else:
        tht2[i] = w[i]*(t[i]-t[i-1]) + tht2[i-1]
    z2 = R2*np.exp(1j*tht2[i])
    sol = root_scalar(calc_tht4,x0=x0,bracket=bracket)
    tht4[i] = sol.root
    z4 = R4*np.exp(1j*tht4[i])
    z5 = z1 + z4
    z3 = z5 - z2
    tht3[i] = np.angle(z3)
    
# positions of A and B
z_A = R2*np.exp(1j*tht2)
z_B = z_A + R3*np.exp(1j*tht3)

plt.figure()
plt.plot(np.real(z_A),np.imag(z_A),label='Point A')
plt.plot(np.real(z_B),np.imag(z_B),label='Point B')
plt.axis('equal')
plt.legend()

# spline curve fits
x_A_spline = interpolate.splrep(t,np.real(z_A))   # spline of point A x-position array
y_A_spline = interpolate.splrep(t,np.imag(z_A))   # spline of point A y-position array
x_B_spline = interpolate.splrep(t,np.real(z_B))   # spline of point B x-position array
y_B_spline = interpolate.splrep(t,np.imag(z_B))   # spline of point B y-position array

# plot position vs. time
_, ax = plt.subplots(2,1)

ax[0].plot(t,interpolate.splev(t,x_A_spline),label='Point A')
ax[0].plot(t,interpolate.splev(t,x_B_spline),label='Point B')
ax[0].set_ylabel('x-position [cm]')
ax[0].set_xlabel('Time [s]')
ax[0].legend()

ax[1].plot(t,interpolate.splev(t,y_A_spline),label='Point A')
ax[1].plot(t,interpolate.splev(t,y_B_spline),label='Point B')
ax[1].set_ylabel('y-position [cm]')
ax[1].set_xlabel('Time [s]')
ax[1].legend()

# plot velocity vs. time
_, ax = plt.subplots(2,1)

ax[0].plot(t,interpolate.splev(t,x_A_spline,der=1),label='Point A')
ax[0].plot(t,interpolate.splev(t,x_B_spline,der=1),label='Point B')
ax[0].set_ylabel('x-velocity [cm/s]')
ax[0].set_xlabel('Time [s]')
ax[0].legend()

ax[1].plot(t,interpolate.splev(t,y_A_spline,der=1),label='Point A')
ax[1].plot(t,interpolate.splev(t,y_B_spline,der=1),label='Point B')
ax[1].set_ylabel('y-velocity [cm/s]')
ax[1].set_xlabel('Time [s]')
ax[1].legend()

# plot acceleration vs. time
_, ax = plt.subplots(2,1)

ax[0].plot(t,interpolate.splev(t,x_A_spline,der=2),label='Point A')
ax[0].plot(t,interpolate.splev(t,x_B_spline,der=2),label='Point B')
ax[0].set_ylabel('x-acceleration [cm/s^2]')
ax[0].set_xlabel('Time [s]')
ax[0].legend()

ax[1].plot(t,interpolate.splev(t,y_A_spline,der=2),label='Point A')
ax[1].plot(t,interpolate.splev(t,y_B_spline,der=2),label='Point B')
ax[1].set_ylabel('y-acceleration [cm/s^2]')
ax[1].set_xlabel('Time [s]')
ax[1].legend()

# compute tangential and normal components of acceleration of point A
V_A = interpolate.splev(t,x_A_spline,der=1) + interpolate.splev(t,y_A_spline,der=1)*1j
A_A = interpolate.splev(t,x_A_spline,der=2) + interpolate.splev(t,y_A_spline,der=2)*1j
a_At = (np.real(V_A)*np.real(A_A)+np.imag(V_A)*np.imag(A_A))/np.abs(V_A)
A_At = a_At*np.exp(1j*np.angle(V_A))
A_An = A_A - A_At
a_An = np.abs(A_An)

plt.figure()
plt.plot(t,a_At,label='Point A tangential acceleration')
plt.plot(t,a_An,label='Point A normal acceleration')
plt.ylabel('Acceleration [cm/s^2]')
plt.xlabel('Time [s]')
plt.legend()

plt.show()