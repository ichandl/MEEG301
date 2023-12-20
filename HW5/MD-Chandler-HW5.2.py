# -*- coding: utf-8 -*-
# Isaac Chandler
# F23 Machine Design HW5.2
# Uses code from "4bar_acceleration.py" and "crank_slider_velocity.py"

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# crank position
w = 10.472                          # rotation rate [rad/s]
t = np.linspace(0,4*np.pi/w,101)     # time array [s]
R2 = 30                              # radius [cm]
z_A = R2*np.exp(1j*w*t)              # complex position [cm]

# slider position
R3 = 100                                   # coupler length [cm]
D = 40                                       # slider vertical offset [cm]
tht3 = np.arcsin((D-R2*np.sin(w*t))/R3)     # coupler angle [rad]
z_B = z_A + R3*np.exp(1j*tht3)              # complex position [cm]

# spline curve fits
x_A_spline = interpolate.splrep(t,np.real(z_A))   # spline of crank x-position array
y_A_spline = interpolate.splrep(t,np.imag(z_A))   # spline of crank y-position array
x_B_spline = interpolate.splrep(t,np.real(z_B))   # spline of slider x-position array
y_B_spline = interpolate.splrep(t,np.imag(z_B))   # spline of slider y-position array

# (a) plot acceleration vs. time
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

# (b) compute tangential and normal components of acceleration of point A
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

# (c) plot acceleration vs. time of point C
z_C = z_A + R3/2*np.exp(1j*tht3)
x_C_spline = interpolate.splrep(t,np.real(z_C))   # spline of slider x-position array
y_C_spline = interpolate.splrep(t,np.imag(z_C))   # spline of slider y-position array

_, ax = plt.subplots(2,1)

ax[0].plot(t,interpolate.splev(t,x_C_spline,der=2),label='Point C')
ax[0].set_ylabel('x-acceleration [cm/s^2]')
ax[0].set_xlabel('Time [s]')
ax[0].legend()

ax[1].plot(t,interpolate.splev(t,y_C_spline,der=2),label='Point C')
ax[1].set_ylabel('y-acceleration [cm/s^2]')
ax[1].set_xlabel('Time [s]')
ax[1].legend()

# (d) Plot Normal, Tangent Acceleration of C
V_C = interpolate.splev(t,x_C_spline,der=1) + interpolate.splev(t,y_C_spline,der=1)*1j
A_C = interpolate.splev(t,x_C_spline,der=2) + interpolate.splev(t,y_C_spline,der=2)*1j
a_Ct = (np.real(V_C)*np.real(A_C)+np.imag(V_C)*np.imag(A_C))/np.abs(V_C)
A_Ct = a_Ct*np.exp(1j*np.angle(V_C))
A_Cn = A_C - A_Ct
a_Cn = np.abs(A_Cn)


Unit = 1*np.exp(1j*tht3)
U_C = np.full(101, Unit)
A_C2 = np.zeros(101)

for i, j in enumerate(A_C):
    A_C2[i] = np.dot(A_C[i],U_C[i])

a_Ct2 = np.real(A_C2)
a_Cn2 = np.imag(A_C2)

plt.figure()
plt.plot(t,a_Ct2,label='Point C tangential acceleration')
plt.plot(t,a_Cn2,label='Point C normal acceleration')
plt.ylabel('Acceleration [cm/s^2]')
plt.xlabel('Time [s]')
plt.legend()

plt.show()