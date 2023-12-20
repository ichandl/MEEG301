# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 19:13:55 2020

@author: Adam Wickenheiser
"""

import numpy as np

def check_Grashof(a,b,c,d):
    # Prints "Grashof","non-Grashof", or "Special case Grashof"
    # based on link lengths a,b,c,d (in any order)
    
    S = min(a,b,c,d)
    L = max(a,b,c,d)
    SL = S + L
    PQ = a + b + c + d - SL
    if SL < PQ:
        print('Grashof')
    elif SL == PQ:
        print('Special case Grashof')
    else:
        print('non-Grashof')
        

# 2-position crank rocker synthesis

# givens
tht4 = 60*np.pi/180
phi = 60*np.pi/180

# check relationship between R2 and R4
R4 = 6
R2 = R4*np.sin(phi/2)
M = R4*(np.exp(1j*(phi+tht4))-np.exp(1j*tht4))
print(2*R2)
print(np.abs(M))

# compute position of crank-ground joint
R3 = 7
K = (R2+R3)/np.abs(M)
z_B1 = R4*np.exp(1j*tht4)
z_O2 = z_B1 + K*M
print('x2 = ',np.real(z_O2))
print('y2 = ',np.imag(z_O2))

# check Grashof condition
R1 = np.abs(z_O2)
check_Grashof(R1, R2, R3, R4)
