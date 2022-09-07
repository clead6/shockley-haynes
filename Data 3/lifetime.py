# -*- coding: utf-8 -*-
"""
Created on Thu Mar 31 10:00:31 2022

@author: h02845cd
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm


def linear(x,m,p):
    return m*x+p

def gaussian(x, a, b, c, d):
    return a*np.exp(-(x - b)**2/(2*c**2)) + d

def chi_squared(observed,expected,error):
    return np.sum((observed-expected)**2/error**2)

#files = ["17V.CSV","19V.CSV","21V.CSV","23V.CSV","25V.CSV","27V.CSV","29V.CSV",
#         "31V.CSV","33V.CSV","35V.CSV","37V.CSV","39V.CSV","41V.CSV","43V.CSV",
#         "45V.CSV","47V.CSV"]
#keys = ["17V","19V","21V","23V","25V","27V","29V","31V","33V","35V","37V",
#        "39V","41V","43V","45V","47V"]
#voltage = np.array([17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47])

files = ["16V.CSV","17V.CSV","18V.CSV",
         "19V.CSV","20V.CSV","21V.CSV",
         "22V.CSV","23V.CSV","24V.CSV","25V.CSV",
         "26V.CSV","27V.CSV","28V.CSV","29V.CSV"]
keys = ["16V","17V","18V",
         "19V","20V","21V",
         "22V","23V","24V","25V",
         "26V","27V","28V","29V"]
voltage = np.array([16,17,18,19,20,21,22,23,24,25,26,27,28,29])


counter=16
data={}
time = np.array([])
voltage = np.array([])
for i in files:
    data["{0}V".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=1
    
max_volt=np.array([])
t_0=np.array([])
unc_t_0=np.array([])
fwhm=np.array([])
unc_fwhm=np.array([])
chi2 = np.array([])
n_0=np.array([])
unc_n_0=np.array([])

count=0 
counter=16
for key in keys:
    
    if count<2:
        time=data[key][1000:2450,3]
        volt=data[key][1000:2450,4]

    else:
        time=data[key][500:2450,3]
        volt=data[key][500:2450,4]
    

    #unc_volt = 0.0005
    unc_volt = volt/np.sqrt(len(volt))

    max_volt=np.append(max_volt,max(volt))
    indexes=np.where(volt==max_volt[count])
    
    if len(indexes[0])==1:
        t_0=np.append(t_0,time[indexes[0]])
        unc_t_0= np.append(unc_t_0,time[1]-time[0])
    
    else:
        time_range = np.array([])
        last_index=indexes[0][-1]
        for i in range(indexes[0][0],last_index-1,1):
            time_range=np.append(time_range, time[i])
        t_0 = np.append(t_0,np.mean(time_range))
        unc_t_0= np.append(unc_t_0,time_range[-1]-time_range[0])

    popt_gau, pcov_gau = curve_fit(gaussian, time, volt, p0=(7, 2e-5, 1e-5, 8))
    n_0=np.append(n_0,popt_gau[0])
    unc_popt_gau = np.sqrt(np.diag(pcov_gau))
    unc_n_0=np.append(unc_n_0, unc_popt_gau[0])
    
    chi2 = np.append(chi2, chi_squared(volt, gaussian(time, popt_gau[0],popt_gau[1],popt_gau[2],popt_gau[3]), unc_volt)/len(volt))

    # fwhm = np.append(fwhm, 2*np.sqrt(2*np.log(2))*popt_gau[2])
    # unc_fwhm = np.append(unc_fwhm, fwhm[count]*unc_popt_gau[2]/popt_gau[2])
    #unc_fwhm = np.append(unc_fwhm, 2*np.sqrt(2*np.log(2))*unc_popt_gau[2])
    #print(fwhm[count], "+/-", unc_fwhm[count])

    #plt.plot(time, gaussian(time, popt_gau[0],popt_gau[1],popt_gau[2],popt_gau[3]), marker=".")
    #plt.scatter(time, volt, marker='.', color="r")
    #plt.xlabel("time")
    #plt.ylabel("voltage")
    #plt.savefig("gaussian{0}".format(counter),dpi=300)
    #plt.show()

    count+=1
    counter+=1
    

x_values = t_0
y_values = np.log(n_0*np.sqrt(t_0))
unc_y = np.sqrt((unc_t_0 / (2*t_0))**2 + (unc_n_0 / n_0)**2)

fit = curve_fit(linear,x_values,y_values)
x_axis = np.linspace(1.7e-5,2.8e-5,500)

chi2_linear = chi_squared(y_values, linear(x_values,fit[0][0],fit[0][1]), unc_y)
red_chi2_linear = chi2_linear / (len(y_values) - 2)

plt.errorbar(x_values, y_values, yerr=unc_y, ls="none")
plt.scatter(x_values,y_values, marker=".")
plt.plot(x_axis,linear(x_axis,fit[0][0],fit[0][1]))
plt.xlabel(r"$t_0 [s]$")
plt.ylabel(r"$\ln({n_0 \sqrt{t_0}})$")
plt.savefig("lifetime", dpi=300)
plt.show()

tau = -1/fit[0][0]