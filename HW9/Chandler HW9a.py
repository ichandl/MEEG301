# Isaac Chandler
# F23 Machine Design
# HW9.1a, uses code from "engine_forces_torques.py"

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

# geometric parameters
w = 1500/60*2*np.pi                  # rotation rate [rad/s]
t = np.linspace(0,2*np.pi/w,1001)    # time array [s]
R2 = .03469                            # crank length [m]
R3 = .11465                            # connecting rod length [m]

# angle splines [rad]
tht2 = w*t                                 # crank angle array [rad]
tht3 = np.arcsin(-R2/R3*np.sin(w*t))       # connecting rod angle [rad]
tht2_spline = interpolate.splrep(t,tht2)   # spline of theta 2 array
tht3_spline = interpolate.splrep(t,tht3)   # spline of theta 3 array

# plot link angles
plt.figure()
plt.plot(t,180/np.pi*tht2,label='tht_2')
plt.plot(t,180/np.pi*tht3,label='tht_3')
plt.ylabel('Angle [deg]')
plt.xlabel('Time [s]')
plt.legend()

# angular accelerations needed for force calculations [rad/s^2]
alpha2 = interpolate.splev(t,tht2_spline,der=2)
alpha3 = interpolate.splev(t,tht3_spline,der=2)

# joint positions
z_A = R2*np.exp(1j*w*t)                  # position of point A [m]
z_B = z_A + R3*np.exp(1j*tht3)           # position of point B [m]

# CG positions
z_CG2 = R2/2*np.exp(1j*tht2)
z_CG3 = R2*np.exp(1j*tht2) + R3/2*np.exp(1j*tht3)
z_CG4 = z_B

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

# displacements
z_12 = -z_CG2
z_32 = z_A - z_CG2
z_23 = z_A - z_CG3
z_43 = z_B - z_CG3
z_34 = z_B - z_CG4

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

# link masses
# rho = 8000             # density of steel [kg/m^3]
# b = 0.010              # link width [m]
# h = 0.002              # link thickness [m]
m2 = .12125#rho*R2*b*h        # mass of crank [kg]
m3 = .24966#rho*R3*b*h        # mass of connecting rod [kg]
m4 = .91864#0.8               # mass of piston [kg]

print('\n')
print('m2 = ',m2,'kg')
print('m3 = ',m3,'kg')
print('m4 = ',m4,'kg')

# link moments of inertia about CG
I_G2 = .00012163#m2/12*(R2**2 + b**2)         # moment of inertia of crank [kg-m^2].00012163#
I_G3 = .00124745#m3/12*(R3**2 + b**2)         # moment of inertia of connecting rod [kg-m^2].00124745#

print('I_G2 = ',I_G2,'kg-m^2')
print('I_G3 = ',I_G3,'kg-m^2')

# applied external force [N]
Fg = np.zeros_like(t)
comb_ratio = 0.25                        # ratio of the cycle to apply Fg
Fg[t < comb_ratio*2*np.pi/w] = 0     # gas force [N]

plt.figure()
plt.plot(t,Fg)
plt.ylabel('Gas Force [N]')
plt.xlabel('Time [s]')

# initialize force/torque arrays
F12x = np.zeros_like(t)
F12y = np.zeros_like(t)
F32x = np.zeros_like(t)
F32y = np.zeros_like(t)
F43x = np.zeros_like(t)
F43y = np.zeros_like(t)
F14y = np.zeros_like(t)
T12 = np.zeros_like(t)

# loop through time array
for i in range(t.size):
    # set up and solve linear system A*X=B
    A = [[1, 0, 1, 0, 0, 0, 0, 0],
         [0, 1, 0, 1, 0, 0, 0, 0],
         [-R12y[i], R12x[i], -R32y[i], R32x[i], 0, 0, 0, 1],
         [0, 0, -1, 0, 1, 0, 0, 0],
         [0, 0, 0, -1, 0, 1, 0, 0],
         [0, 0, R23y[i], -R23x[i], -R43y[i], R43x[i], 0, 0],
         [0, 0, 0, 0, -1, 0, 0, 0],
         [0, 0, 0, 0, 0, -1, 1, 0]]
    B = [m2*aCG2x[i], m2*aCG2y[i], I_G2*alpha2[i],
         m3*aCG3x[i], m3*aCG3y[i], I_G3*alpha3[i],
         m4*aCG4x[i] + Fg[i], m4*aCG4y[i]]
    X = np.linalg.solve(A, B)

    # assign solution to force/torque arrays
    F12x[i] = X[0]
    F12y[i] = X[1]
    F32x[i] = X[2]
    F32y[i] = X[3]
    F43x[i] = X[4]
    F43y[i] = X[5]
    F14y[i] = X[6]
    T12[i] = X[7]

# plot results
_, ax = plt.subplots(2,1)

ax[0].plot(t,F12x,label='F12')
ax[0].plot(t,F32x,label='F32')
ax[0].plot(t,F43x,label='F43')
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