# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 18:24:43 2021

@author: Adam Wickenheiser
"""
# Modified by Isaac Chandler
# F23 MEEG301 HW7

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('/Users/esackchan/Desktop/Code/Junior Fall/Machine Design/HW7/MD_HW7B_Torque.csv',header=1)
t = np.array(data.iloc[:,0])
T = np.array(data.iloc[:,1])
# F12x = np.array(data.iloc[:,1])
# F23x = np.array(data.iloc[:,2])
# F34x = np.array(data.iloc[:,3])
# F41x = np.array(data.iloc[:,4])



plt.plot(t,T)
plt.legend()
plt.title("Problem 2 Torque; Solidworks")
plt.xlabel("Time (t)")
plt.ylabel("Torque (N*mm)")

plt.show()