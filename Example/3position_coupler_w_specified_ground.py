# -*- coding: utf-8 -*-
"""
Created on Tue Sep 22 12:49:17 2020

@author: Adam Wickenheiser
"""

import numpy as np
from scipy.optimize import shgo

def rigid_body_position_update(z_C1,z_D1,z_C2,z_D2,z_E1,z_F1):
    # Given absolute complex positions of C1,D1,C2,D2,E1,F1, compute E2,F2
    
    # compute translation (not needed) and rotation from position 1 to 2
    # z21 = z_C2 - z_C1
    z_D1C1 = z_D1 - z_C1
    z_D2C2 = z_D2 - z_C2
    tht21 = np.angle(z_D2C2) - np.angle(z_D1C1)
    
    # relative positions of E1 and F1 with respect to C1
    z_E1C1 = z_E1 - z_C1
    z_F1C1 = z_F1 - z_C1
    
    # rotate these relative positions through angle tht2
    z_E2C2 = z_E1C1*np.exp(1j*tht21)
    z_F2C2 = z_F1C1*np.exp(1j*tht21)
    
    # compute absolute positions of E2 and F2
    z_E2 = z_C2 + z_E2C2
    z_F2 = z_C2 + z_F2C2
    
    return z_E2, z_F2
        
def intersect_perp_bisectors(z_A1,z_B1,z_A2,z_B2):
    # Find the intersection of the perpendicular bisectors between A1 and A2 and between 
    # B1 and B2
    
    # find midpoints
    z_Am = 0.5*(z_A1+z_A2)
    z_Bm = 0.5*(z_B1+z_B2)
    
    # compute angles of bisectors
    thtA = np.angle(z_A2-z_A1) + np.pi/2
    thtB = np.angle(z_B2-z_B1) + np.pi/2
    
    # solve a*x = b
    a = [[np.cos(thtA), -np.cos(thtB)],[np.sin(thtA), -np.sin(thtB)]]
    b = [np.real(z_Bm-z_Am), np.imag(z_Bm-z_Am)]
    x = np.linalg.solve(a, b)
    
    z_O = z_Am + x[0]*np.exp(1j*thtA)
    
    return z_O
    
def dist_to_ground_joints(x):
    # Computes the location of the rotopole using input 
    # x = [real(z_E1C1), imag(z_E1C1), real(z_F1C1), imag(z_F1C1)],
    # then calculates its distance to the desired locations z_O2,z_O4.
    
    z_E1C1 = x[0] + x[1]*1j
    z_F1C1 = x[2] + x[3]*1j
    z_E1 = z_C1 + z_E1C1
    z_F1 = z_C1 + z_F1C1
    
    z_E2, z_F2 = rigid_body_position_update(z_C1,z_D1,z_C2,z_D2,z_E1,z_F1)
    z_E3, z_F3 = rigid_body_position_update(z_C1,z_D1,z_C3,z_D3,z_E1,z_F1)
    
    z_O2_guess = intersect_perp_bisectors(z_E1,z_E2,z_E2,z_E3)
    z_O4_guess = intersect_perp_bisectors(z_F1,z_F2,z_F2,z_F3)
    
    d = np.abs(z_O2_guess - z_O2) + np.abs(z_O4_guess - z_O4)
    return d
    
    


# 3-position coupler output synthesis with specified ground joints

# givens
z_C1 = 2.164 + 1.26j
z_D1 = z_C1 + np.exp(1j*30*np.pi/180)
z_C2 = z_C1 - 1.236 + 2.138j
z_D2 = z_C2 + np.exp(1j*-32.5*np.pi/180)
z_C3 = z_C1 - 2.5 + 2.931j
z_D3 = z_C3 + np.exp(1j*-69.8*np.pi/180)
z_O2 = 0 + 0j
z_O4 = 2.164 + 2.19 + 0j
print('z_C1 = ',z_C1)
print('z_D1 = ',z_D1)
print('z_C2 = ',z_C2)
print('z_D2 = ',z_D2)
print('z_C3 = ',z_C3)
print('z_D3 = ',z_D3)
print('z_O2 = ',z_O2)
print('z_O4 = ',z_O4)

# solve the optimization problem
res = shgo(dist_to_ground_joints, bounds=((-5,5),(-5,5),(-5,5),(-5,5)),
           options={'f_min': 0.001, 'maxev': 100})
print(res.message)
print(res.fun)

z_E1C1 = res.x[0] + res.x[1]*1j
z_F1C1 = res.x[2] + res.x[3]*1j
z_E1 = z_C1 + z_E1C1
z_F1 = z_C1 + z_F1C1
print('z_E1 = ',z_E1)
print('z_F1 = ',z_F1)


