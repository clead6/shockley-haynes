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

log = np.log(np.sqrt(t_0)*abs(n_0))


popt_mo, pcov_mo = curve_fit(proportional, voltage/l, d/t_0)
unc_popt_mo = np.sqrt(np.diag(pcov_mo))
mobility = popt_mo[0]
unc_mobility = unc_popt_mo[0]
print("Mobility = ", mobility, "+/-", unc_mobility, "cm^2 s^-1 V-1, rel = ", unc_mobility/mobility)

x_data = np.linspace(0,25,1000)
y_unc = (d/t_0)*np.sqrt((unc_d/d)**2+(unc_t_0/t_0)**2)
plt.plot(x_data, proportional(x_data, popt_mo[0]))
plt.scatter(voltage/l, d/t_0, marker=".")
plt.errorbar(voltage/l, d/t_0, yerr=y_unc,ls='none')
plt.xlabel(r"$V_s/l$ [V cm$^{-1}$]")
plt.ylabel(r"$d/t$ [cm s$^{-1}$]")
plt.title("Mobility for d=0.3cm varying $V_s$")
plt.savefig("mobility_fpb",dpi=300)
plt.show()


popt_dif, pcov_dif = curve_fit(proportional, t_0**3, fwhm**2)
unc_popt_dif = np.sqrt(np.diag(pcov_dif))
diffusivity = popt_dif[0]*0.3**2/(16*np.log(2))
unc_diffusivity = diffusivity*np.sqrt((unc_popt_dif[0]/popt_dif[0])**2+(2*unc_d/d)**2)
print("Diffusivity = ", diffusivity, "+/-", unc_diffusivity, "cm^2 s^-1, rel = ", unc_diffusivity/diffusivity)

x_data = np.linspace(0, 1.9e-14, 1000)
y_unc = 2*fwhm*unc_fwhm
plt.plot(x_data, proportional(x_data, popt_dif[0]))
plt.scatter(t_0**3, fwhm**2, marker=".")
plt.errorbar(t_0**3, fwhm**2, yerr=y_unc, ls="none")
plt.xlabel(r"$t_0^3$ [s$^3$]")
plt.ylabel(r"$t_p^2$ [s$^3$]")
plt.title("Diffusivity for d=0.3cm varying $V_s$")
plt.savefig("diffusion_fpb",dpi=300)
plt.show()


popt_lif, pcov_lif = curve_fit(linear, t_0, np.log(abs(n_0)*np.sqrt(t_0)))
unc_popt_lif = np.sqrt(np.diag(pcov_lif))
lifetime = -1/popt_lif[0]
unc_lifetime = lifetime*abs(unc_popt_lif[0]/popt_lif[0])
print("Lifetime = ", lifetime, "+/-", unc_lifetime, "s, rel = ", unc_lifetime/lifetime)

x_data = np.linspace(1.2e-5,2.6e-5,1000)
y_unc = np.sqrt((unc_t_0/(2*t_0))**2+(unc_n_0/n_0)**2)
plt.plot(x_data, linear(x_data, popt_lif[0], popt_lif[1]))
plt.scatter(t_0, np.log(abs(n_0)*np.sqrt(t_0)), marker=".")
plt.errorbar(t_0, np.log(abs(n_0)*np.sqrt(t_0)), yerr=y_unc, ls="none")
plt.xlabel(r"$t_0$ [s]")
plt.ylabel(r"$\ln({n_0 \sqrt{t_0}})$")
plt.title("Lifetime for d=0.3cm varying $V_s$")
plt.savefig("lifetime_fpb", dpi=300)
plt.show()

# mobility correction analysis
mobility_array = (d*l) / (t_0*voltage)
aveg_mob = np.average(mobility_array)
x = (2*constants.k*300*l) / (constants.e*voltage*d) * (t_0 / lifetime + 0.5)
mobility_correction = mobility_array * (np.sqrt(1 + x**2) - x)
plt.scatter(voltage, mobility_array, marker=".", label="measured values")
plt.scatter(voltage, mobility_correction, marker=".", label="corrected values")
plt.legend()
plt.ylim([0,1450])
plt.xlabel(r"$V_s$ [V]")
plt.ylabel(r"$\mu$ [cm$^2$ s$^{-1}$ V$^{-1}$]")
plt.title(r"Mobility for d=0.3cm varying $V_s$")
plt.savefig("mobility_horizontal_fpb", dpi=300)
plt.show()

mobility_corrected = np.average(mobility_correction)
mobility_corrected_unc = np.std(mobility_correction)
print("Mobility final = ", mobility_corrected, "+/-", mobility_corrected_unc, "cm^2 s^-1 V-1, rel = ", mobility_corrected_unc/mobility_corrected)

# einstein ratio
einstein_ratio = diffusivity / mobility
print("Einstein relation = ", einstein_ratio, "V")
