# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:46:58 2022

@author: clead
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy import constants

def proportional(x, m):
    return m*x

def linear(x, m, p):
    return m*x + p

def gaussian(x, a, b, c, d):
    return a*np.exp(-(x - b)**2/(2*c**2)) + d

def exponential(x, b, c, d):
    return np.exp(-b*x + c) + d

def chi_squared(observed,expected,error):
    return np.sum((observed-expected)**2/error**2)

def non_linear(x,a,b):
    return a*np.exp(x*b)/np.sqrt(x)

files = ["17V.CSV","18V.CSV",
         "19V.CSV","20V.CSV","21V.CSV",
         "22V.CSV","23V.CSV","24V.CSV","25V.CSV",
         "26V.CSV","27V.CSV","28V.CSV","29V.CSV",
         "30V.CSV","31V.CSV","32V.CSV","33V.CSV",
         "34V.CSV","35V.CSV","36V.CSV","37V.CSV",
         "38V.CSV","39V.CSV","40V.CSV","41V.CSV",
         "42V.CSV","43V.CSV","44V.CSV","45V.CSV",
         "46V.CSV","47V.CSV","48V.CSV",
         "49V.CSV","50V.CSV"]
keys = ["17V","18V",
         "19V","20V","21V",
         "22V","23V","24V","25V",
         "26V","27V","28V","29V",
         "30V","31V","32V","33V",
         "34V","35V","36V","37V",
         "38V","39V","40V","41V",
         "42V","43V","44V","45V",
         "46V","47V","48V",
         "49V","50V"]
voltage = np.array([17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
                    40,41,42,43,44,45,46,47,48,49,50])

counter=17
data={}
for i in files:
    data["{0}V".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=1


max_volt=np.array([])
t_0=np.array([])
unc_t_0=np.array([])
fwhm=np.array([])
unc_fwhm=np.array([])
n_0=np.array([])
unc_n_0=np.array([])
chi2 = np.array([])
count=0 
counter=17
for key in keys:

    if 17<=counter<=22:

        time = data[key][400:2400,3]
        volt = data[key][400:2400,4] 
        peaks, properties = find_peaks(volt, distance=1500, width=200, height = -5)
    
        volt_left = volt[:int(properties["left_ips"])-250]
        volt_right = volt[int(properties["right_ips"])+250:]
        volt_concat = np.concatenate((volt_left,volt_right))
        time_concat = time[:len(volt_concat)]
    
    if 22<counter<35:
        time = data[key][300:2000,3]
        volt = data[key][300:2000,4]

        peaks, properties = find_peaks(volt, distance=1500, width=200, height = -5)
    
        volt_left = volt[:int(properties["left_ips"])-290]
        volt_right = volt[int(properties["right_ips"])+200:]
        volt_concat = np.concatenate((volt_left,volt_right))
        time_concat = time[:len(volt_concat)]        
        
    if counter>34:
        time = data[key][300:1700,3]
        volt = data[key][300:1700,4]        

        peaks, properties = find_peaks(volt, distance=1500, width=200, height = -5)
        
        volt_left = volt[:int(properties["left_ips"])-250]
        volt_right = volt[int(properties["right_ips"])+250:]
        volt_concat = np.concatenate((volt_left,volt_right))
        time_concat = time[:len(volt_concat)]
        
        
    popt_exp, pcov_exp = curve_fit(exponential, time_concat, volt_concat, p0=(10, -1, -1.5),check_finite=False, maxfev=10000000)
    background = exponential(time, popt_exp[0],popt_exp[1],popt_exp[2])
    volt_subt = volt - (background-background[-1])    


    peaks, properties = find_peaks(volt_subt, distance=1500, width=200, height = -5)     

    widths = properties["widths"]
    t_0 = np.append(t_0, time[peaks])
    n_0 = np.append(n_0, 2*(properties["width_heights"]-properties["peak_heights"]))        
    fwhm = np.append(fwhm,(time[int(properties["right_ips"])] - time[int(properties["left_ips"])]))
    
    unc_t_0 = np.append(unc_t_0, time[peaks]*0.01)
    unc_n_0 = np.append(unc_n_0, 2*(properties["width_heights"]-properties["peak_heights"])*0.01)
    unc_fwhm = np.append(unc_fwhm, 0.01*(time[int(properties["right_ips"])] - time[int(properties["left_ips"])]))
    
    
    # plt.scatter(time, volt_subt-properties["peak_heights"]+properties["prominences"], marker='.', color="r")
    # plt.plot(time[peaks], volt[peaks]-properties["peak_heights"]+properties["prominences"], "x")
    # plt.vlines(x=time[peaks], ymin=volt[peaks] -properties["peak_heights"], ymax = volt[peaks]-properties["peak_heights"]+properties["prominences"], color = "C1")
    # plt.hlines(y=properties["width_heights"]-properties["peak_heights"]+properties["prominences"], xmin=time[int(properties["left_ips"])], xmax=time[int(properties["right_ips"])], color = "C1")
    # plt.xlabel("time [s]")
    # plt.ylabel("voltage [V]")
    # plt.title("Signal for d=0.3cm and V$_s$={0}V".format(counter))
    # plt.savefig("gaussian{0}_fpb".format(counter),dpi=300)
    # plt.show()

    count+=1
    counter+=1
 
d=0.3 #in cm
unc_d=0.005 #in cm
l=2.0 #in cm
unc_l=0.1 #in cm

unc_voltage=0.1 # in V


popt_nlin, pcov_nlin = curve_fit(non_linear, t_0, abs(n_0), p0=(0.02,-100000))
unc_popt_nlin = np.sqrt(np.diag(pcov_nlin))
lifetime = -1/popt_nlin[1]
unc_lifetime = lifetime*abs(unc_popt_nlin[1]/popt_nlin[1])
print("Lifetime = ", lifetime, "+/-", unc_lifetime, "s, rel = ", unc_lifetime/lifetime)

x_data = np.linspace(1.2e-5,2.6e-5,1000)
plt.scatter(t_0,abs(n_0),marker='.')
plt.plot(x_data, non_linear(x_data, popt_nlin[0], popt_nlin[1]))
plt.errorbar(t_0, abs(n_0), yerr=unc_n_0,ls="none")
plt.xlabel(r"$t_0$ [s]")
plt.ylabel(r"$n_0$")
plt.title("Lifetime for d=0.3cm varying $V_s$")
plt.savefig("lifetime_non_linear", dpi=300)
plt.show()
