# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 10:12:49 2022

@author: h02845cd
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def linear(x, m):
    return m*x

def gaussian(x, a, b, c, d):
    return a*np.exp(-(x - b)**2/(2*c**2)) + d

def exponential(x, a, b, c, d):
    return d*np.exp(-a*x + b) + c

def full_fit(x, a, b, c, d, e, f, g, h):
    return gaussian(x, a, b, c, d) + exponential(x, e,f,g,h)

def chi_squared(observed,expected,error):
    return np.sum((observed-expected)**2/error**2)

files = ["17V.CSV","18V.CSV",
         "19V.CSV","20V.CSV","21V.CSV",
         "22V.CSV","23V.CSV","24V.CSV","25V.CSV",
         "26V.CSV","27V.CSV","28V.CSV","29V.CSV"]
keys = ["17V","18V",
         "19V","20V","21V",
         "22V","23V","24V","25V",
         "26V","27V","28V","29V"]
voltage = np.array([17,18,19,20,21,22,23,24,25,26,27,28,29])


counter=17
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
count=0 
counter=17
for key in keys:
    if (counter<25):
        time = data[key][400:2450,3]
        volt = data[key][400:2450,4]
    else:
        time = data[key][500:2000,3]
        volt = data[key][500:2000,4]

    unc_volt = volt/np.sqrt(len(volt))

    max_volt=np.append(max_volt,max(volt))
    indexes=np.where(volt==max_volt[count])
    
    if len(indexes[0])==1:
        t_0=np.append(t_0,time[indexes[0]])
        # unc_t_0= np.append(unc_t_0,time[1]-time[0])
    
    else:
        time_range = np.array([])
        last_index=indexes[0][-1]
        for i in range(indexes[0][0],last_index-1,1):
            time_range=np.append(time_range, time[i])
        t_0 = np.append(t_0,np.mean(time_range))
        # unc_t_0= np.append(unc_t_0,time_range[-1]-time_range[0])
        
    if counter<25:
    
        min_index = np.where(volt==min(volt[:int(len(volt)/2)]))[0][0]
        popt_exp, pcov_exp = curve_fit(exponential, time[:min_index], volt[:min_index], p0=(10, -1, -1.5, 1),maxfev=100000)
        #popt_exp, pcov_exp = curve_fit(exponential, time, volt, p0=(10, -1, -1.5, 1),maxfev=1000000)
    
        # popt_full, pcov_full = curve_fit(full_fit, time, volt, p0=(7, 2e-5, 1e-5, 8, 10, -1, -10, 1),maxfev=100000)
        # plt.plot(time, full_fit(time, popt_full[0],popt_full[1],popt_full[2],popt_full[3], 
        #                         popt_full[4],popt_full[5],popt_full[6],popt_full[7]), marker=".")
    
        #popt_exp_2, pcov_exp_2 = curve_fit(exponential, time[2000:], volt[2000:], p0=(10, -1, -1.5, 1),maxfev=10000)
    
        plt.plot(time, exponential(time, popt_exp[0],popt_exp[1],popt_exp[2],popt_exp[3]), marker=".")
        plt.scatter(time, volt, marker='.', color="r")
        plt.xlabel("time")
        plt.ylabel("voltage")
        plt.title("{0}V".format(counter))
        plt.savefig("gaussian{0}_background".format(counter),dpi=300)
        plt.show()
        
    
        background = exponential(time, popt_exp[0],popt_exp[1],popt_exp[2],popt_exp[3])
        #background_2 = exponential(time, popt_exp_2[0],popt_exp_2[1],popt_exp_2[2],popt_exp_2[3])
        volt_subt = volt - (background-background[-1]) 
        # volt_subt = volt - (background-background[-1]) - (background_2-background_2[-1])
        popt_gau, pcov_gau = curve_fit(gaussian, time, volt_subt, p0=(7, 2e-5, 1e-5, 8))
        plt.scatter(time, volt_subt, marker='.', color="r")
        
    else:
        popt_gau, pcov_gau = curve_fit(gaussian, time, volt, p0=(7, 2e-5, 1e-5, 8))
        plt.scatter(time, volt, marker='.', color="r")

    plt.plot(time, gaussian(time, popt_gau[0], popt_gau[1], popt_gau[2], popt_gau[3]), marker=".")

    plt.xlabel("time")
    plt.ylabel("voltage")
    plt.title("{0}V".format(counter))
    plt.savefig("gaussian{0}".format(counter),dpi=300)
    plt.show()
    
    count+=1
    counter+=1
    
    unc_popt_gau = np.sqrt(np.diag(pcov_gau))


# popt, pcov = curve_fit(linear,t_0**3, fwhm**2)
# unc_popt = np.sqrt(np.diag(pcov))

# x_data = np.linspace(0, 2.1e-14, 1000)
# plt.plot(x_data, linear(popt[0], x_data))

# plt.scatter(t_0**3, fwhm**2, marker=".")
# #plt.errorbar(t_0**3, fwhm**2, yerr=2*fwhm**2*unc_popt_gau[2]/popt_gau[2], ls="none")
# plt.xlabel("t_0^3")
# plt.ylabel("t_p^2")
# plt.savefig("diffusion constant - day 3",dpi=300)
# plt.show()


# d=0.3 #in cm
# unc_d=0.01 #in cm

# dif = popt[0]*0.3**2/(16*np.log(2))
# unc_dif=dif*np.sqrt((unc_popt[0]/popt[0])**2+(2*unc_d/d)**2)
# print("diffusion constant=", dif, "+/-", unc_dif, ", ", unc_dif/dif*100, "%")

