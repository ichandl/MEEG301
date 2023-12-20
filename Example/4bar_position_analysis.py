# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 16:54:36 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

link_lengths = [8, 5, 7, 7]           # link lengths [cm]
tht2 = 70*np.pi/180                     # angle of crank [rad]
z1 = link_lengths[0]                    # ground link
z2 = link_lengths[1]*np.exp(1j*tht2)    # crank link
link_config = 'open'                    # 'open' or 'crossed'

# choose which solution of tht4 to look for
tht_div = np.angle(z2-z1)
if link_config == 'open':
    if tht2 >= 0:
        bracket = [tht_div-np.pi,tht_div]      # angle range to search within
        x0 = tht_div-np.pi/2                   # initial guess for tht4
    else:
        bracket = [tht_div,tht_div+np.pi]      # angle range to search within
        x0 = tht_div+np.pi/2                   # initial guess for tht4
else:  # 'crossed'
    if tht2 >= 0:
        bracket = [tht_div,tht_div+np.pi]      # angle range to search within
        x0 = tht_div+np.pi/2                   # initial guess for tht4
    else:
        bracket = [tht_div-np.pi,tht_div]      # angle range to search within
        x0 = tht_div-np.pi/2                   # initial guess for tht4

# function whose root we want to find within the bracket
def calc_tht4(tht4_guess):
    z4 = link_lengths[3]*np.exp(1j*tht4_guess)
    z5 = z1 + z4
    z3 = z5 - z2
    r = np.abs(z3)
    return r - link_lengths[2]    # r - R3

# find the root
sol = root_scalar(calc_tht4,x0=x0,bracket=bracket)  # calls the calc_tht4 function repeatedly with new guesses
tht4 = sol.root
z4 = link_lengths[3]*np.exp(1j*tht4)
z5 = z1 + z4
z3 = z5 - z2
tht3 = np.angle(z3)
print('theta 3 = ',tht3*180/np.pi,'deg.')
print('theta 4 = ',tht4*180/np.pi,'deg.')

# plot the resulting linkage configuration
linkage = [0, z2, z5, z1, 0]
plt.plot(np.real(linkage),np.imag(linkage),lw=2,marker='.',ms=16)
plt.axis('equal')
plt.show();