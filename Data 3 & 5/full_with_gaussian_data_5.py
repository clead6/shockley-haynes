# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 11:46:58 2022

@author: clead
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

def exponential(x, b, c, d):
    return np.exp(-b*x + c) + d

def chi_squared(observed,expected,error):
    return np.sum((observed-expected)**2/error**2)

def non_linear(x,a,b):
    return a*np.exp(x*b)/np.sqrt(x)

files = ["3.0mm.CSV","3.1mm.CSV",
         "3.2mm.CSV","3.3mm.CSV","3.4mm.CSV",
         "3.5mm.CSV","3.6mm.CSV","3.7mm.CSV","3.8mm.CSV",
         "3.9mm.CSV","4.0mm.CSV","4.1mm.CSV","4.2mm.CSV",
         "4.3mm.CSV","4.4mm.CSV","4.5mm.CSV","4.6mm.CSV",
         "4.7mm.CSV","4.8mm.CSV","4.9mm.CSV","5.0mm.CSV",
         "5.1mm.CSV","5.2mm.CSV","5.3mm.CSV","5.4mm.CSV",
         "5.5mm.CSV"]
keys = ["30","31","32","33","34",
         "35","36","37","38",
         "39","40","41","42",
         "43","44","45","46",
         "47","48","49","50",
         "51","52","53","54",
         "55"]
d = np.array([3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5])
unc_d=0.05

counter=30
data={}
for i in files:
    data["{0}".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
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
counter=30
for key in keys:

    if 29<counter<36:

        time = data[key][800:1700,3]
        volt = data[key][800:1700,4] 
        
        peaks, properties = find_peaks(volt, distance=1500, width=100, height = -8)
    
        volt_left = volt[:int(properties["left_ips"])-150]
        volt_right = volt[int(properties["right_ips"])+150:]
        volt_middle = np.full(int(properties["right_ips"])+150-(int(properties["left_ips"])-150), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time
    
    if 36<=counter<45:
        time = data[key][850:1800,3]
        volt = data[key][850:1800,4] 
        
        peaks, properties = find_peaks(volt, distance=1000, width=100, height = -8)
    
        volt_left = volt[:int(properties["left_ips"])-150]
        volt_right = volt[int(properties["right_ips"])+150:]
        volt_middle = np.full(int(properties["right_ips"])+150-(int(properties["left_ips"])-150), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time
        
    if 45<=counter<50:
        time = data[key][900:2000,3]
        volt = data[key][900:2000,4] 
        
        peaks, properties = find_peaks(volt, distance=1000, width=50, height = -8)
    
        volt_left = volt[:int(properties["left_ips"])-150]
        volt_right = volt[int(properties["right_ips"])+150:]
        volt_middle = np.full(int(properties["right_ips"])+150-(int(properties["left_ips"])-150), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time
        
    if 50<=counter:
        time = data[key][1000:2000,3]
        volt = data[key][1000:2000,4] 
        
        peaks, properties = find_peaks(volt, distance=800, width=50, height = -8)
    
        volt_left = volt[:int(properties["left_ips"])-150]
        volt_right = volt[int(properties["right_ips"])+150:]
        volt_middle = np.full(int(properties["right_ips"])+150-(int(properties["left_ips"])-150), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time
        
    popt_exp, pcov_exp = curve_fit(exponential, time_concat, volt_concat, p0=(10, -1, -1.5),check_finite=False, maxfev=10000000)
    background = exponential(time, popt_exp[0],popt_exp[1],popt_exp[2])
    volt_subt = volt - (background-background[-1])    
    
    popt_gau, pcov_gau = curve_fit(gaussian, time, volt_subt,  p0=(7, 2e-5, 1e-5, 8))
    unc_popt_gau = np.sqrt(np.diag(pcov_gau))

    t_0 = np.append(t_0, popt_gau[1])
    n_0 = np.append(n_0, popt_gau[0])        
    fwhm = np.append(fwhm, popt_gau[2]*2*np.sqrt(2*np.log(2)))
    
    unc_t_0 = np.append(unc_t_0, unc_popt_gau[1])
    unc_n_0 = np.append(unc_n_0, unc_popt_gau[0])
    unc_fwhm = np.append(unc_fwhm, unc_popt_gau[2]*2*np.sqrt(2*np.log(2)))
    
    plt.scatter(time, volt_subt-popt_gau[3], marker='.', color="r")
    plt.plot(time, gaussian(time, popt_gau[0], popt_gau[1], popt_gau[2], popt_gau[3])-popt_gau[3])
    plt.xlabel("time [s]")
    plt.ylabel("voltage [V]")
    plt.title("Signal for V$_s$=30V and d={:.1f}mm".format(counter*0.1))
    plt.savefig("gaussian{0}".format(counter),dpi=300)
    plt.show()

    count+=1
    counter+=1
 
l=2.0 #in cm
unc_l=0.1 #in cm
voltage = 30.0
unc_voltage = 0.1

d*=0.1 #in cm
unc_d*=0.1 # in cm

log = np.log(np.sqrt(t_0)*abs(n_0))

popt_mo, pcov_mo = curve_fit(proportional, voltage/d, l/t_0)
unc_popt_mo = np.sqrt(np.diag(pcov_mo))
mobility = popt_mo[0]
unc_mobility = unc_popt_mo[0]
print("Mobility = ", mobility, "+/-", unc_mobility, "cm^2 s^-1 V-1, rel = ", unc_mobility/mobility)

x_data = np.linspace(0,100,1000)
y_unc = (l/t_0)*np.sqrt((unc_l/l)**2+(unc_t_0/t_0)**2)
plt.plot(x_data, proportional(x_data, popt_mo[0]))
plt.scatter(voltage/d, l/t_0, marker=".")
plt.errorbar(voltage/d, l/t_0, yerr=y_unc,ls='none')
plt.xlabel(r"$V_s/d$ [V cm$^{-1}$]")
plt.ylabel(r"$l/t$ [cm s$^{-1}$]")
plt.title("Mobility for V$_s$=30V varying d")
plt.savefig("mobility_fpb",dpi=300)
plt.show()


popt_dif, pcov_dif = curve_fit(proportional, t_0**3, (d*fwhm)**2)
unc_popt_dif = np.sqrt(np.diag(pcov_dif))
diffusivity = popt_dif[0]/(16*np.log(2))
unc_diffusivity = diffusivity*(unc_popt_dif[0]/popt_dif[0])
print("Diffusivity = ", diffusivity, "+/-", unc_diffusivity, "cm^2 s^-1, rel = ", unc_diffusivity/diffusivity)

x_data = np.linspace(0, 3.0e-14, 1000)
y_unc = 2*(d*fwhm)**2*np.sqrt((unc_fwhm/fwhm)**2+(unc_d/d)**2)
plt.plot(x_data, proportional(x_data, popt_dif[0]))
plt.scatter(t_0**3, (d*fwhm)**2, marker=".")
plt.errorbar(t_0**3, (d*fwhm)**2, yerr=y_unc, ls="none")
plt.xlabel(r"$t_0^3$ [s$^3$]")
plt.ylabel(r"$(t_p d)^2$ [s$^3$ cm$^2$]")
plt.title("Diffusivity for V$_s$=30V varying d")
plt.savefig("diffusion_fpb",dpi=300)
plt.show()


popt_lif, pcov_lif = curve_fit(linear, t_0, np.log(np.sqrt(t_0)*abs(n_0)))
unc_popt_lif = np.sqrt(np.diag(pcov_lif))
lifetime = -1/popt_lif[0]
unc_lifetime = lifetime*abs(unc_popt_lif[0]/popt_lif[0])
print("Lifetime = ", lifetime, "+/-", unc_lifetime, "s, rel = ", unc_lifetime/lifetime)

x_data = np.linspace(1.7e-5,3.1e-5,1000)
y_unc = np.sqrt((unc_t_0/(2*t_0))**2+(unc_n_0/n_0)**2)
plt.plot(x_data, linear(x_data, popt_lif[0], popt_lif[1]))
plt.scatter(t_0, np.log(np.sqrt(t_0)*abs(n_0)), marker=".")
plt.errorbar(t_0, np.log(np.sqrt(t_0)*abs(n_0)), yerr=y_unc, ls="none")
plt.xlabel(r"$t_0$ [s]")
plt.ylabel(r"$\ln({n_0 \sqrt{t_0}})$")
plt.title("Lifetime for d=0.3cm varying $V_s$")
plt.savefig("lifetime_fpb", dpi=300)
plt.show()

# non linear lifetime
popt_nlin, pcov_nlin = curve_fit(non_linear, t_0, abs(n_0), p0=(0.02,-100000))
unc_popt_nlin = np.sqrt(np.diag(pcov_nlin))
lifetime = -1/popt_nlin[1]
unc_lifetime = lifetime*abs(unc_popt_nlin[1]/popt_nlin[1])
print("Lifetime (non linear) = ", lifetime, "+/-", unc_lifetime, "s, rel = ", unc_lifetime/lifetime)

x_data = np.linspace(1.7e-5,3.3e-5,1000)
plt.scatter(t_0,abs(n_0),marker='.')
plt.plot(x_data, non_linear(x_data, popt_nlin[0], popt_nlin[1]))
plt.errorbar(t_0, abs(n_0), yerr=unc_n_0,ls="none")
plt.xlabel(r"$t_0$ [s]")
plt.ylabel(r"$n_0$")
plt.title(r"Lifetime for $V_s$=30 V varying $d$")
plt.savefig("lifetime_non_linear", dpi=300)
plt.show()

# mobility vs distance
mobility_array = (d*l) / (t_0*voltage)
mobility_array_unc = mobility_array*np.sqrt((unc_d/d)**2+(unc_l/l)**2+(unc_t_0/t_0)**2+(unc_voltage/voltage)**2)
plt.scatter(d, mobility_array, marker=".")
plt.errorbar(d,mobility_array,yerr=mobility_array_unc, ls='none')
plt.ylim([0,1400])
plt.xlabel(r"$d$ [cm]")
plt.ylabel(r"$\mu$ [cm$^2$ s$^{-1}$ V$^{-1}$]")
plt.title(r"Mobility for V$_s$=30V varying d")
plt.savefig("mobility_horizontal_fpb", dpi=300)
plt.show()

# eintein ratio
einstein_ratio = diffusivity / mobility
print("Einstein relation = ", einstein_ratio, "V")
