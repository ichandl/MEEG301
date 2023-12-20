# -*- coding: utf-8 -*-
# Isaac Chandler
# F23 Machine Design HW6.1
# Uses code from "4bar_forces_single_config.py"

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
tht0 = 60*np.pi/180         # instantaneous crank angle [rad]
delta = -31*np.pi/180         # angle between AP and AB [rad]
w0 = 10                      # instantaneous crank velocity [rad/s]
alpha0 = 5                  # instantaneous crank acceleration [rad/s^2]
rocker_dir = 'up'            # rocker direction 'up' or 'down'
z1 = R1                      # ground link

# set up arrays for simulation
t = np.linspace(0,2,501)                                     # time array [s]
tht2 = alpha0/2*t**2 + (w0-alpha0)*t + tht0 - w0 + alpha0/2  # crank angle array [rad]
tht3 = np.zeros_like(t)                                      # coupler angle array [rad]
tht4 = np.zeros_like(t)                                      # rocker angle array [rad]

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

# angles, angular velocities, and angular accelerations needed for force calculations
theta2 = interpolate.splev(1.0,tht2_spline)
omega2 = interpolate.splev(1.0,tht2_spline,der=1)
alpha2 = interpolate.splev(1.0,tht2_spline,der=2)

theta3 = interpolate.splev(1.0,tht3_spline)
omega3 = interpolate.splev(1.0,tht3_spline,der=1)
alpha3 = interpolate.splev(1.0,tht3_spline,der=2)

theta4 = interpolate.splev(1.0,tht4_spline)
omega4 = interpolate.splev(1.0,tht4_spline,der=1)
alpha4 = interpolate.splev(1.0,tht4_spline,der=2)

print('theta2:',180/np.pi*theta2,'deg.')
print('omega2:',omega2,'rad/s')
print('alpha2:',alpha2,'rad/s^2')

print('theta3:',180/np.pi*theta3,'deg.')
print('omega3:',omega3,'rad/s')
print('alpha3:',alpha3,'rad/s^2')

print('theta4:',180/np.pi*theta4,'deg.')
print('omega4:',omega4,'rad/s')
print('alpha4:',alpha4,'rad/s^2')

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

# CG accelerations needed for force calculations [m/s^2]
aCG2x = interpolate.splev(1.0,x_CG2_spline,der=2)
aCG2y = interpolate.splev(1.0,y_CG2_spline,der=2)
aCG3x = interpolate.splev(1.0,x_CG3_spline,der=2)
aCG3y = interpolate.splev(1.0,y_CG3_spline,der=2)
aCG4x = interpolate.splev(1.0,x_CG4_spline,der=2)
aCG4y = interpolate.splev(1.0,y_CG4_spline,der=2)

print('\n')
print('aCG2x:',aCG2x,'m/s')
print('aCG2y:',aCG2y,'m/s')
print('aCG3x:',aCG3x,'m/s')
print('aCG3y:',aCG3y,'m/s')
print('aCG4x:',aCG4x,'m/s')
print('aCG4y:',aCG4y,'m/s')

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
R12x = interpolate.splev(1.0,R12x_spline)
R12y = interpolate.splev(1.0,R12y_spline)
R32x = interpolate.splev(1.0,R32x_spline)
R32y = interpolate.splev(1.0,R32y_spline)
R23x = interpolate.splev(1.0,R23x_spline)
R23y = interpolate.splev(1.0,R23y_spline)
R43x = interpolate.splev(1.0,R43x_spline)
R43y = interpolate.splev(1.0,R43y_spline)
R34x = interpolate.splev(1.0,R34x_spline)
R34y = interpolate.splev(1.0,R34y_spline)
R14x = interpolate.splev(1.0,R14x_spline)
R14y = interpolate.splev(1.0,R14y_spline)
RP3x = interpolate.splev(1.0,RP3x_spline)
RP3y = interpolate.splev(1.0,RP3y_spline)

print('\n')
print('R12x:',R12x,'m')
print('R12y:',R12y,'m')
print('R32x:',R32x,'m')
print('R32y:',R32y,'m')
print('R23x:',R23x,'m')
print('R23y:',R23y,'m')
print('R43x:',R43x,'m')
print('R43y:',R43y,'m')
print('R34x:',R34x,'m')
print('R34y:',R34y,'m')
print('R14x:',R14x,'m')
print('R14y:',R14y,'m')
print('RP3x:',RP3x,'m')
print('RP3y:',RP3y,'m')

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

# form matrices and solve A*X = B
FPx = 0          # x-force at P [N]
FPy = -100       # y-force at P [N]
T4 = 0           # torque at O4 [N-m]
A = [[1, 0, 1, 0, 0, 0, 0, 0, 0],
     [0, 1, 0, 1, 0, 0, 0, 0, 0],
     [-R12y, R12x, -R32y, R32x, 0, 0, 0, 0, 1],
     [0, 0, -1, 0, 1, 0, 0, 0, 0],
     [0, 0, 0, -1, 0, 1, 0, 0, 0],
     [0, 0, R23y, -R23x, -R43y, R43x, 0, 0, 0],
     [0, 0, 0, 0, -1, 0, 1, 0, 0],
     [0, 0, 0, 0, 0, -1, 0, 1, 0],
     [0, 0, 0, 0, R34y, -R34x, -R14y, R14x, 0]]
B = [m2*aCG2x, m2*aCG2y, I_G2*alpha2, m3*aCG3x-FPx, m3*aCG3y-FPy, I_G3*alpha3-RP3x*FPy+RP3y*FPx,
     m4*aCG4x, m4*aCG4y, I_G4*alpha4-T4]
X = np.linalg.solve(A, B)

# convert to polar coordinates
F12 = X[0] + X[1]*1j
F32 = X[2] + X[3]*1j
F43 = X[4] + X[5]*1j
F14 = X[6] + X[7]*1j
T12 = X[8]

print('\n')
print('F12 = ',np.abs(F12),'N @',np.angle(F12)*180/np.pi,'deg.')
print('F32 = ',np.abs(F32),'N @',np.angle(F32)*180/np.pi,'deg.')
print('F43 = ',np.abs(F43),'N @',np.angle(F43)*180/np.pi,'deg.')
print('F14 = ',np.abs(F14),'N @',np.angle(F14)*180/np.pi,'deg.')
print('T12 = ',T12,'N-m')