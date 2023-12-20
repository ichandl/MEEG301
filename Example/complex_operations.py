# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 17:28:51 2020

@author: Adam Wickenheiser
"""

import numpy as np

# rotation by 30 deg.
z1 = 2+1j
z2 = z1*np.exp(1j*30*np.pi/180)
print('\nRotation example')
print(z2)
print('z1 angle = ',np.angle(z1)*180/np.pi,'deg.')
print('z2 angle = ',np.angle(z2)*180/np.pi,'deg.')

# four-bar linkage, all angles known
z2 = 8*np.exp(1j*40*np.pi/180)
z3 = 6*np.exp(1j*-10.82*np.pi/180)
z5 = z2 + z3
print('\nFour-bar linkage, all angles known example')
print('x2 = ',np.real(z2))
print('y2 = ',np.imag(z2))
print('x5 = ',np.real(z5))
print('y5 = ',np.imag(z5))
print(15+5*np.exp(1j*126.56*np.pi/180))

# four-bar linkage, all coordinates known
z2 = 6.36 + 6.36j
R2 = np.abs(z2)
tht2 = np.angle(z2)*180/np.pi
z5 = 13.15 + 4.65j
z3 = z5 - z2
R3 = np.abs(z3)
tht3 = -np.angle(z3)*180/np.pi
z1 = 15
z4 = z5 - z1
R4 = np.abs(z4)
tht4 = np.angle(z4)*180/np.pi
print('\nFour-bar linkage, all coordinates known example')
print('R2 = ',R2)
print('theta2 = ',tht2)
print('z3 = ',z3)
print('R3 = ',R3)
print('theta3 = ',tht3)
print('z4 = ',z4)
print('R4 = ',R4)
print('theta4 = ',tht4)