# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:42:30 2022

@author: h02845cd
"""

import numpy as np
import matplotlib.pyplot as plt

voltage = np.genfromtxt("diode_x.csv", dtype=float, delimiter=",")[:,4]
current = np.genfromtxt("diode_y.csv", dtype=float, delimiter=",")[:,4]

plt.scatter(voltage, current, marker=".")
plt.xlabel("Voltage [V]")
plt.ylabel("Current")
plt.title("Diode-behaviour of the Ge Bar")
plt.savefig("Diode-behaviour", dpi=300)
plt.show()
