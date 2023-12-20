# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 13:42:08 2020

@author: Adam Wickenheiser
"""

import numpy as np

def rigid_body_position_update(z_C1,z_D1,z_C2,z_D2,z_E1,z_F1):
    # Given absolute complex positions of C1,D1,C2,D2,E1,F1, compute E2,F2
    
    # compute translation (not needed) and rotation from position 1 to 2
    # z21 = z_C2 - z_C1
    z_D1C1 = z_D1 - z_C1
    z_D2C2 = z_D2 - z_C2
    tht21 = np.angle(z_D2C2) - np.angle(z_D1C1)
    
    # relative posotions of E1 and F1 with respect to C1
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
    
    


# 3-position coupler output synthesis

z_C1 = -1.5 + 9.58j
z_D1 = 1.89 + 11.70j
z_C2 = 2.72 + 10.31j
z_D2 = 6.55 + 9.19j
z_C3 = 5.96 + 8.00j
z_D3 = 9.45 + 9.96j

z_E1 = -0.04 + 9.24j
z_F1 = 1.84 + 8.02j

z_E2, z_F2 = rigid_body_position_update(z_C1,z_D1,z_C2,z_D2,z_E1,z_F1)
print('z_E2 = ',z_E2)
print('z_F2 = ',z_F2)
z_E3, z_F3 = rigid_body_position_update(z_C1,z_D1,z_C3,z_D3,z_E1,z_F1)
print('z_E3 = ',z_E3)
print('z_F3 = ',z_F3)

z_O2 = intersect_perp_bisectors(z_E1,z_E2,z_E2,z_E3)
print('z_O2 = ',z_O2)
z_O4 = intersect_perp_bisectors(z_F1,z_F2,z_F2,z_F3)
print('z_O4 = ',z_O4)
