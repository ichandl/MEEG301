# -*- coding: utf-8 -*-
# Isaac Chandler
# F23 Machine Design HW6.2
# Uses code from "4bar_forces_full_cycle.py"

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.optimize import root_scalar
from scipy import interpolate

# givens
R1 = 2.22                     # ground link length [m]
R2 = 1.0                    # crank length [m]
R3 = 2.06                     # coupler length [m]
R4 = 2.33                     # rocker length [m]
AP = 3.06                     # distance between points A and P [m]
delta = -31*np.pi/180         # angle between AP and AB [rad]
w = 15                        # crank velocity [rad/s]
n = 1                        # number of cycles to plot
rocker_dir = 'up'            # rocker direction 'up' or 'down'
z1 = R1                      # ground link

# set up arrays for simulation
t = np.linspace(0,n*2*np.pi/w,501)   # time array [s]
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


# loop through all time steps and calculate angles
for i in range(t.size):
    z2 = R2*np.exp(1j*tht2[i])
    sol = root_scalar(calc_tht4,x0=x0,bracket=bracket)
    tht4[i] = sol.root
    z4 = R4*np.exp(1j*tht4[i])
    z5 = z1 + z4
    z3 = z5 - z2
    tht3[i] = np.angle(z3)

# angle splines [rad]
tht2 = np.unwrap(tht2)
tht3 = np.unwrap(tht3)
tht4 = np.unwrap(tht4)
tht2_spline = interpolate.splrep(t,tht2)   # spline of theta 2 array
tht3_spline = interpolate.splrep(t,tht3)   # spline of theta 3 array
tht4_spline = interpolate.splrep(t,tht4)   # spline of theta 4 array

# plot link angles
plt.figure()
plt.plot(t,180/np.pi*tht2,label='tht_2')
plt.plot(t,180/np.pi*tht3,label='tht_3')
plt.plot(t,180/np.pi*tht4,label='tht_4')
plt.ylabel('Angle [deg]')
plt.xlabel('Time [s]')
plt.legend()

# angular accelerations needed for force calculations [rad/s^2]
alpha2 = interpolate.splev(t,tht2_spline,der=2)
alpha3 = interpolate.splev(t,tht3_spline,der=2)
alpha4 = interpolate.splev(t,tht4_spline,der=2)

# CG positions
z_CG2 = R2/2*np.exp(1j*tht2)
z_CG3 = R2*np.exp(1j*tht2) + R3/2*np.exp(1j*tht3)
z_CG4 = R1 + R4/2*np.exp(1j*tht4)

plt.figure()
plt.plot(np.real(z_CG2),np.imag(z_CG2),label='Link 2')
plt.plot(np.real(z_CG3),np.imag(z_CG3),label='Link 3')
plt.plot(np.real(z_CG4),np.imag(z_CG4),label='Link 4')
plt.title('CG positions')
plt.axis('equal')
plt.legend()

# CG position splines
x_CG2_spline = interpolate.splrep(t,np.real(z_CG2))   # spline of CG_2 x-position array
y_CG2_spline = interpolate.splrep(t,np.imag(z_CG2))   # spline of CG_2 y-position array
x_CG3_spline = interpolate.splrep(t,np.real(z_CG3))   # spline of CG_3 x-position array
y_CG3_spline = interpolate.splrep(t,np.imag(z_CG3))   # spline of CG_3 y-position array
x_CG4_spline = interpolate.splrep(t,np.real(z_CG4))   # spline of CG_4 x-position array
y_CG4_spline = interpolate.splrep(t,np.imag(z_CG4))   # spline of CG_4 y-position array

# CG accelerations needed for force calculations [m/s^2]
aCG2x = interpolate.splev(t,x_CG2_spline,der=2)
aCG2y = interpolate.splev(t,y_CG2_spline,der=2)
aCG3x = interpolate.splev(t,x_CG3_spline,der=2)
aCG3y = interpolate.splev(t,y_CG3_spline,der=2)
aCG4x = interpolate.splev(t,x_CG4_spline,der=2)
aCG4y = interpolate.splev(t,y_CG4_spline,der=2)

# joint positions
z_A = R2*np.exp(1j*tht2)
z_B = z_A + R3*np.exp(1j*tht3)
z_P = z_A + AP*np.exp(1j*tht3)

# displacements
z_12 = -z_CG2
z_32 = z_A - z_CG2
z_23 = z_A - z_CG3
z_43 = z_B - z_CG3
z_34 = z_B - z_CG4
z_14 = z1 - z_CG4
z_P3 = z_P - z_CG3

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
RP3x_spline = interpolate.splrep(t,np.real(z_P3))
RP3y_spline = interpolate.splrep(t,np.imag(z_P3))

# displacements needed for force calculations [m]
R12x = interpolate.splev(t,R12x_spline)
R12y = interpolate.splev(t,R12y_spline)
R32x = interpolate.splev(t,R32x_spline)
R32y = interpolate.splev(t,R32y_spline)
R23x = interpolate.splev(t,R23x_spline)
R23y = interpolate.splev(t,R23y_spline)
R43x = interpolate.splev(t,R43x_spline)
R43y = interpolate.splev(t,R43y_spline)
R34x = interpolate.splev(t,R34x_spline)
R34y = interpolate.splev(t,R34y_spline)
R14x = interpolate.splev(t,R14x_spline)
R14y = interpolate.splev(t,R14y_spline)
RP3x = interpolate.splev(t,RP3x_spline)
RP3y = interpolate.splev(t,RP3y_spline)

# link masses
rho_steel = 7800            # density of steel [kg/m^3]
rho_aluminum = 2800         # density of aluminum [kh/m^3]
b = 0.050              # link width [m]
h = 0.025              # link thickness [m]
m2 = rho_steel*R2*b*h        # mass of link 2 [kg]
m3 = rho_steel*R3*b*h        # mass of link 3 [kg]
m4 = rho_steel*R4*b*h        # mass of link 4 [kg]
m_plate = 0.5*rho_aluminum*R3*AP*math.sin(delta)*b       #mass of the aluminum plate

print('\n')
print('m2 = ',m2,'kg')
print('m3 = ',m3,'kg')
print('m4 = ',m4,'kg')
print('m_plate = ',m_plate,'kg')

# link moments of inertia about CG [kg-m^2]
I_G2 = m2/12*(R2**2 + b**2)
I_G3 = m3/12*(R3**2 + b**2)
I_G4 = m4/12*(R4**2 + b**2)
I_plate = (m_plate/6)*(R3**2+(AP*math.sin(delta))**2)

print('I_G2 = ',I_G2,'kg-m^2')
print('I_G3 = ',I_G3,'kg-m^2')
print('I_G4 = ',I_G4,'kg-m^2')
print('I_plate = ',I_plate,'kg-m^2')

# applied external force [N]
FP = 20*np.exp(1j*-np.pi/2)*np.ones_like(t)
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
         m3*aCG3x[i]-FPx[i], m3*aCG3y[i]-FPy[i],
         I_G3*alpha3[i]-RP3x[i]*FPy[i]+RP3y[i]*FPx[i],
         m4*aCG4x[i], m4*aCG4y[i], I_G4*alpha4[i]]
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
ax[0].set_ylabel('x-force [N]')
ax[0].set_xlabel('Time [s]')
ax[0].legend()

ax[1].plot(t,F12y,label='F12')
ax[1].plot(t,F32y,label='F32')
ax[1].plot(t,F43y,label='F43')
ax[1].plot(t,F14y,label='F14')
ax[1].set_ylabel('y-force [N]')
ax[1].set_xlabel('Time [s]')
ax[1].legend()

plt.figure()
plt.plot(t,T12,label='T12')
plt.ylabel('Torque [N-m]')
plt.xlabel('Time [s]')
plt.show()