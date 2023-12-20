# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 19:37:17 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar
from scipy import interpolate

# givens
R1 = 9.625                           # ground link length [in]
R2 = 2.0                             # crank length [in]
R3 = 8.375                           # coupler length [in]
R4 = 7.187                           # rocker length [in]
BP = 3.75                            # distance between points B and P [in]
w = 500/60*2*np.pi                   # crank velocity [rad/s]
rocker_dir = 'up'                    # rocker direction 'up' or 'down'
z1 = R1                              # ground link

# set up arrays for simulation
t = np.linspace(0,2*np.pi/w,501)     # time array [s]
tht2 = w*t                           # crank angle array [rad]
tht3 = np.zeros_like(t)              # coupler angle array [rad]
tht4 = np.zeros_like(t)              # rocker angle array [rad]

# choose which solution of tht4 to look for
if rocker_dir == 'up':
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
    z2 = R2*np.exp(1j*tht2[i])
    sol = root_scalar(calc_tht4,x0=x0,bracket=bracket)
    tht4[i] = sol.root
    z4 = R4*np.exp(1j*tht4[i])
    z5 = z1 + z4
    z3 = z5 - z2
    tht3[i] = np.angle(z3)

# angle splines
tht2_spline = interpolate.splrep(t,tht2)   # spline of theta 2 array
tht3_spline = interpolate.splrep(t,tht3)   # spline of theta 3 array
tht4_spline = interpolate.splrep(t,tht4)   # spline of theta 4 array

# angular accelerations needed for force calculations
alpha2 = interpolate.splev(t,tht2_spline,der=2)
alpha3 = interpolate.splev(t,tht3_spline,der=2)
alpha4 = interpolate.splev(t,tht4_spline,der=2)

# CG positions
z_CG2 = R2/2*np.exp(1j*tht2)
z_CG3 = R2*np.exp(1j*tht2) + R3/2*np.exp(1j*tht3)
z_CG4 = R1 + R4/2*np.exp(1j*tht4)

# CG position splines
x_CG2_spline = interpolate.splrep(t,np.real(z_CG2))   # spline of CG_2 x-position array
y_CG2_spline = interpolate.splrep(t,np.imag(z_CG2))   # spline of CG_2 y-position array
x_CG3_spline = interpolate.splrep(t,np.real(z_CG3))   # spline of CG_3 x-position array
y_CG3_spline = interpolate.splrep(t,np.imag(z_CG3))   # spline of CG_3 y-position array
x_CG4_spline = interpolate.splrep(t,np.real(z_CG4))   # spline of CG_4 x-position array
y_CG4_spline = interpolate.splrep(t,np.imag(z_CG4))   # spline of CG_4 y-position array

# CG accelerations needed for force calculations [ft/s^2]
aCG2x = interpolate.splev(t,x_CG2_spline,der=2)/12
aCG2y = interpolate.splev(t,y_CG2_spline,der=2)/12
aCG3x = interpolate.splev(t,x_CG3_spline,der=2)/12
aCG3y = interpolate.splev(t,y_CG3_spline,der=2)/12
aCG4x = interpolate.splev(t,x_CG4_spline,der=2)/12
aCG4y = interpolate.splev(t,y_CG4_spline,der=2)/12

# link masses
gamma = 0.28           # specific weight of steel [lb/in^3]
g = 32.2               # acceleration due to gravity [ft/s^2]
b = 2                  # link width [in]
h = 1                  # link thickness [in]
m2 = gamma*R2*b*h/g    # mass of link 2 [slug]
m3 = gamma*R3*b*h/g    # mass of link 3 [slug]
m4 = gamma*R4*b*h/g    # mass of link 4 [slug]
mB = 14.5/g            # lumped mass at point B [slug]

print('\n')
print('m2 = ',m2,'slug')
print('m3 = ',m3,'slug')
print('m4 = ',m4,'slug')

# link moments of inertia
I_G2 = m2/12*(R2**2 + b**2)/144
I_G3 = m3/12*(R3**2 + b**2)/144
I_G4 = m4/12*(R4**2 + b**2)/144

print('I_G2 = ',I_G2,'slug-ft^2')
print('I_G3 = ',I_G3,'slug-ft^2')
print('I_G4 = ',I_G4,'slug-ft^2')

# joint positions
z_A = R2*np.exp(1j*tht2)
z_B = z_A + R3*np.exp(1j*tht3)
z_P = R1 + (R4+BP)*np.exp(1j*tht4)

# joint position splines
x_B_spline = interpolate.splrep(t,np.real(z_B))   # spline of point B x-position array
y_B_spline = interpolate.splrep(t,np.imag(z_B))   # spline of point B y-position array

# joint accelerations needed for force calculations [ft/s^2]
aBx = interpolate.splev(t,x_B_spline,der=2)/12
aBy = interpolate.splev(t,y_B_spline,der=2)/12

# displacements
z_12 = -z_CG2
z_32 = z_A - z_CG2
z_23 = z_A - z_CG3
z_43 = z_B - z_CG3
z_34 = z_B - z_CG4
z_14 = z1 - z_CG4
z_P4 = z_P - z_CG4

# displacement splines
R12x_spline = interpolate.splrep(t,np.real(z_12))
R12y_spline = interpolate.splrep(t,np.imag(z_12))
R32x_spline = interpolate.splrep(t,np.real(z_32))
R32y_spline = interpolate.splrep(t,np.imag(z_32))
R23x_spline = interpolate.splrep(t,np.real(z_23))
R23y_spline = interpolate.splrep(t,np.imag(z_23))
R43x_spline = interpolate.splrep(t,np.real(z_43))
R43y_spline = interpolate.splrep(t,np.imag(z_43))
R34x_spline = interpolate.splrep(t,np.real(z_34))
R34y_spline = interpolate.splrep(t,np.imag(z_34))
R14x_spline = interpolate.splrep(t,np.real(z_14))
R14y_spline = interpolate.splrep(t,np.imag(z_14))
RP4x_spline = interpolate.splrep(t,np.real(z_P4))
RP4y_spline = interpolate.splrep(t,np.imag(z_P4))

# displacements needed for force calculations [ft]
R12x = interpolate.splev(t,R12x_spline)/12
R12y = interpolate.splev(t,R12y_spline)/12
R32x = interpolate.splev(t,R32x_spline)/12
R32y = interpolate.splev(t,R32y_spline)/12
R23x = interpolate.splev(t,R23x_spline)/12
R23y = interpolate.splev(t,R23y_spline)/12
R43x = interpolate.splev(t,R43x_spline)/12
R43y = interpolate.splev(t,R43y_spline)/12
R34x = interpolate.splev(t,R34x_spline)/12
R34y = interpolate.splev(t,R34y_spline)/12
R14x = interpolate.splev(t,R14x_spline)/12
R14y = interpolate.splev(t,R14y_spline)/12
RP4x = interpolate.splev(t,RP4x_spline)/12
RP4y = interpolate.splev(t,RP4y_spline)/12

# applied external force
FP = 540*np.exp(1j*(tht4+np.pi/2))
FPx = np.real(FP)
FPy = np.imag(FP)

# initialize force/torque arrays
F12x = np.zeros_like(t)
F12y = np.zeros_like(t)
F32x = np.zeros_like(t)
F32y = np.zeros_like(t)
F43x = np.zeros_like(t)
F43y = np.zeros_like(t)
F14x = np.zeros_like(t)
F14y = np.zeros_like(t)
T12 = np.zeros_like(t)

# loop through time array
for i in range(t.size):
    # set up and solve linear system A*X=B
    A = [[1, 0, 1, 0, 0, 0, 0, 0, 0],
         [0, 1, 0, 1, 0, 0, 0, 0, 0],
         [-R12y[i], R12x[i], -R32y[i], R32x[i], 0, 0, 0, 0, 1],
         [0, 0, -1, 0, 1, 0, 0, 0, 0],
         [0, 0, 0, -1, 0, 1, 0, 0, 0],
         [0, 0, R23y[i], -R23x[i], -R43y[i], R43x[i], 0, 0, 0],
         [0, 0, 0, 0, -1, 0, 1, 0, 0],
         [0, 0, 0, 0, 0, -1, 0, 1, 0],
         [0, 0, 0, 0, R34y[i], -R34x[i], -R14y[i], R14x[i], 0]]
    B = [m2*aCG2x[i], m2*aCG2y[i], I_G2*alpha2[i],
         m3*aCG3x[i], m3*aCG3y[i], I_G3*alpha3[i],
         m4*aCG4x[i]+mB*aBx[i]-FPx[i], m4*aCG4y[i]+mB*aBy[i]-FPy[i],
         I_G4*alpha4[i]-RP4x[i]*FPy[i]+RP4y[i]*FPx[i]+R34x[i]*mB*aBy[i]-R34y[i]*mB*aBx[i]]
    X = np.linalg.solve(A, B)

    # assign solution to force/torque arrays
    F12x[i] = X[0]
    F12y[i] = X[1]
    F32x[i] = X[2]
    F32y[i] = X[3]
    F43x[i] = X[4]
    F43y[i] = X[5]
    F14x[i] = X[6]
    F14y[i] = X[7]
    T12[i] = X[8]

# plot results
_, ax = plt.subplots(2,1)

ax[0].plot(t,F12x,label='F12')
ax[0].plot(t,F32x,label='F32')
ax[0].plot(t,F43x,label='F43')
ax[0].plot(t,F14x,label='F14')
ax[0].set_ylabel('x-force [lb]')
ax[0].set_xlabel('Time [s]')
ax[0].legend()

ax[1].plot(t,F12y,label='F12')
ax[1].plot(t,F32y,label='F32')
ax[1].plot(t,F43y,label='F43')
ax[1].plot(t,F14y,label='F14')
ax[1].set_ylabel('y-force [lb]')
ax[1].set_xlabel('Time [s]')
ax[1].legend()

plt.figure()
plt.plot(t,T12,label='T12')
plt.ylabel('Torque [lb-ft]')
plt.xlabel('Time [s]')