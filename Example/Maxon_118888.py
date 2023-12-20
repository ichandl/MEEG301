# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 17:17:25 2021

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt

# motor constants
V_nom = 18           # nominal voltage [V]
w_nl = 13100         # no-load speed [rpm]
i_nl = 0.349         # no-load current [A]
T_nom = 45.9         # nominal torque [mNm]
w_nom = 11500        # nominal speed [rpm]
T_stall = 407        # stall torque [mNm]
kt = 13              # torque constant [mNm/A]
w_max = 25000        # max speed [rpm]

# loading condition (to be driven through gearhead)
T_load = 500         # load torque [mNm]
w_load = 120         # load speed [rpm]
eta_gear = 0.75      # gearhead efficiency

# contour plot of efficiency
T = np.linspace(0,T_stall,101)      # array of torques to plot over [mNm]
w = np.linspace(0,w_max,101)        # array of speeds to plot over [rpm]
V = np.zeros((T.size,w.size))
eta = np.zeros((T.size,w.size))
Tf = kt*i_nl                        # friction torque [mNm]
for i in range(w.size):
    V[:,i] = V_nom*(T/T_stall + w[i]/w_nl)+1e-6             # calculate voltage at this T and w [V]
    eta[:,i] = T/1000*w[i]*(2*np.pi/60)/(V[:,i]*(T+Tf)/kt)  # calculate efficiency at this T and w

#plt.figure(figsize=(15,15))
plt.figure(figsize=(15,15))
c_V = plt.contour(T, w, np.transpose(V))
c_eta = plt.contour(T, w, np.transpose(eta), linestyles='dashed')
plt.plot([T_nom, T_nom], [0, w_max],'r')

plt.clabel(c_V, c_V.levels, inline=True, fmt='%.2f V', fontsize=18)
plt.clabel(c_eta, c_eta.levels, inline=True, fmt='%.2f', fontsize=18)
plt.title('Maxon 118888 motor', fontsize=28)
plt.xlabel('Torque [mNm]', fontsize=24)
plt.ylabel('Speed [rpm]', fontsize=24)
plt.axis([0, T_stall, 0, w_max])
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)

# plot reflected load line
gear_ratio = np.linspace(1,200,201)
gear_ratio_markers = np.geomspace(1,200,11)
T_input = T_load/gear_ratio/eta_gear
w_input = w_load*gear_ratio
T_input_markers = T_load/gear_ratio_markers/eta_gear
w_input_markers = w_load*gear_ratio_markers
plt.plot(T_input, w_input, color='orange', linewidth=2)
plt.plot(T_input_markers, w_input_markers, linestyle='none', marker='o', mfc='orange', mec='k')
for i in range(T_input_markers.size):
    plt.annotate(f'{gear_ratio_markers[i]:.2f}', (T_input_markers[i],w_input_markers[i]),
                 fontsize=18)

# output coordinates in upper right corner of figure window
def format_coord(x, y):
    V = V_nom*(x/T_stall + y/w_nl)
    eta = x/1000*y*(2*np.pi/60)/(V*(x+Tf)/kt)
    return f'torque = {x:.2f} mNm, speed = {y:.2f} rpm, voltage = {V:.2f} V, efficiency = {eta:.2f}'

plt.gca().format_coord = format_coord

