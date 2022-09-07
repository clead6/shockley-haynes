# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 15:09:30 2022

@author: h02845cd
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def proportional(x, m):
    return m*x

def linear(x, m, p):
    return m*x + p

def gaussian(x, a, b, c, d):
    return a*np.exp(-(x - b)**2/(2*c**2)) + d

def chi_squared(observed,expected,error):
    return np.sum((observed-expected)**2/error**2)

files = ["set1.CSV","set2.CSV","set3.CSV","set4.CSV","set5.CSV"]
keys = ["set1","set2","set3","set4","set5"]
voltage = np.array([17,17,17,17,17])

counter=1
data={}
for i in files:
    data["set{0}".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
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
counter=1
for key in keys:

    time = data[key][400:2450,3]
    volt = data[key][400:2450,4]
        
    peaks, properties = find_peaks(volt, distance=1500, width=200, height = -3)
    widths = properties["widths"]
    t_0 = np.append(t_0, time[peaks])
    n_0 = np.append(n_0, 2*properties["width_heights"])
    fwhm = np.append(fwhm,2*np.sqrt(2*np.log(2))*(time[int(properties["right_ips"])] - time[int(properties["left_ips"])]))
    
    plt.scatter(time, volt, marker='.', color="r")
    plt.plot(time[peaks], volt[peaks], "x")
    plt.vlines(x=time[peaks], ymin=volt[peaks] - properties["prominences"], ymax = volt[peaks], color = "C1")
    plt.hlines(y=properties["width_heights"], xmin=time[int(properties["left_ips"])], xmax=time[int(properties["right_ips"])], color = "C1")
    plt.xlabel("time")
    plt.ylabel("voltage")
    plt.title("17V Set {0}".format(counter))
    plt.savefig("gaussian_17V_set{0}_fp".format(counter),dpi=300)
    plt.show()

    count+=1
    counter+=1
 
d=0.3 #in cm
unc_d=0.005 #in cm
l=2.0 #in cm
unc_l=0.1 #in cm

average_t_0 = np.mean(t_0)
error_t_0 = np.std(t_0)
print("t_0: ", average_t_0, "+/-", error_t_0, "rel=", error_t_0/average_t_0)

average_fwhm = np.mean(fwhm)
error_fwhm = np.std(fwhm)
print("fwhm: ", average_fwhm, "+/-", error_fwhm, "rel=", error_fwhm/average_fwhm)

average_n_0 = np.mean(n_0)
error_n_0 = np.std(n_0)
print("n_0: ", average_n_0, "+/-", error_n_0, "rel=", error_n_0/average_n_0)


# popt_mo, pcov_mo = curve_fit(proportional, voltage/l, d/t_0)
# unc_popt_mo = np.sqrt(np.diag(pcov_mo))
# mobility = popt_mo[0]
# unc_mobility = unc_popt_mo[0]
# print("Mobility = ", mobility, "+/-", unc_mobility, "cm^2 s^-1 V-1")

# x_data = np.linspace(0,15,1000)
# plt.plot(x_data, proportional(x_data, popt_mo[0]))
# plt.scatter(voltage/l, d/t_0, marker=".")
# # plt.errorbar(x_axis,y_axis, yerr=y_unc,ls='none')
# plt.xlabel(r"$V_s/l$")
# plt.ylabel(r"$d/t$")
# plt.title("Mobility")
# plt.savefig("mobility_fp",dpi=300)
# plt.show()


# popt_dif, pcov_dif = curve_fit(proportional, t_0**3, fwhm**2)
# unc_popt_dif = np.sqrt(np.diag(pcov_dif))
# diffusivity = popt_dif[0]*0.3**2/(16*np.log(2))
# unc_diffusivity = diffusivity*np.sqrt((unc_popt_dif[0]/popt_dif[0])**2+(2*unc_d/d)**2)
# print("Diffusivity = ", diffusivity, "+/-", unc_diffusivity, "cm^2 s^-1")

# x_data = np.linspace(0, 2.1e-14, 1000)
# plt.plot(x_data, proportional(x_data, popt_dif[0]))
# plt.scatter(t_0**3, fwhm**2, marker=".")
# # plt.errorbar(t_0**3, fwhm**2, yerr=2*fwhm**2*unc_popt_gau[2]/popt_gau[2], ls="none")
# plt.xlabel(r"$t_0^3 [s]$")
# plt.ylabel(r"$t_p^2 [s]$")
# plt.title("Diffusivity")
# plt.savefig("diffusion_fp",dpi=300)
# plt.show()


# popt_lif, pcov_lif = curve_fit(linear, t_0, np.log(abs(n_0)*np.sqrt(t_0)))
# unc_popt_lif = np.sqrt(np.diag(pcov_lif))
# lifetime = -1/popt_lif[0]
# unc_lifetime = lifetime*abs(unc_popt_lif[0]/popt_lif[0])
# print("Lifetime = ", lifetime, "+/-", unc_lifetime, "s")

# x_data = np.linspace(1.7e-5,2.8e-5,500)
# plt.plot(x_data, linear(x_data, popt_lif[0], popt_lif[1]))
# plt.scatter(t_0, np.log(abs(n_0)*np.sqrt(t_0)), marker=".")
# # plt.errorbar(x_values, y_values, yerr=unc_y, ls="none")
# plt.xlabel(r"$t_0 [s]$")
# plt.ylabel(r"$\ln({n_0 \sqrt{t_0}})$")
# plt.title("Lifetime")
# plt.savefig("lifetime_fp", dpi=300)
# plt.show()
