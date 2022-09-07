# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:42:30 2022

@author: h02845cd
"""

import numpy as np
import matplotlib.pyplot as plt

voltage = np.genfromtxt("resistance_x.csv", dtype=float, delimiter=",")[:,4]
current = np.genfromtxt("resistance_y.csv", dtype=float, delimiter=",")[:,4]

plt.scatter(voltage, current, marker=".")
plt.xlabel("Voltage [V]")
plt.ylabel("Current")
plt.title("Resistance-behaviour of the Ge Bar")
plt.savefig("Resistance-behaviour", dpi=300)
plt.show()
