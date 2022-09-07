# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def linear(m,x):
    return m*x

def chi_squared(observed,expected,error):
    return np.sum((observed-expected)**2/error**2)

voltage = np.array([45,40,35,30,25,20,15])
t_0 = np.array([0.00001394, 0.00001503, 0.00001613, 0.00001766, 0.00001966, 0.0000231,0.0000294])
unc_t_0 = np.array([2E-8,6.2E-07,5.80E-07,9.2E-07,1.440E-06,0.00000212,0.0000034])

rel_unc=unc_t_0/t_0
print(rel_unc)


x_axis = voltage/(20*0.1)
y_axis = 1.79*0.1/t_0
y_unc = y_axis*np.sqrt((unc_t_0/t_0)**2+(0.05/1.79)**2)
fit = curve_fit(linear, x_axis, y_axis,p0=(1000,))

x_plot = np.linspace(0,23,100)

mobility=fit[0][0]
unc_mobility = np.sqrt(np.diag(fit[1]))
print("mobility= ", mobility, "+/-", unc_mobility)


chi2=chi_squared(y_axis, linear(fit[0][0],x_axis), y_unc)
red_chi2=chi2/(len(x_axis)-1)

print(red_chi2)

plt.scatter(x_axis, y_axis,marker=".")
plt.errorbar(x_axis,y_axis, yerr=y_unc,ls='none')
plt.plot(x_plot, fit[0][0]*x_plot)
plt.xlabel("V_s/l")
plt.ylabel("d/t")
plt.savefig("mobility_v1",dpi=300)
plt.show()



