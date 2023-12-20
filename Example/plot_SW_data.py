# -*- coding: utf-8 -*-
"""
Created on Tue Nov  2 18:24:43 2021

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('Forces_x.csv',header=1)
t = np.array(data.iloc[:,0])
F12x = np.array(data.iloc[:,1])
F23x = np.array(data.iloc[:,2])
F34x = np.array(data.iloc[:,3])
F41x = np.array(data.iloc[:,4])


plt.plot(t,F12x,t,F23x,t,F34x,t,F41x)