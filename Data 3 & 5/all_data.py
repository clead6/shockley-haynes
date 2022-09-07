# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:00:44 2022

@author: h02845cd
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

files_volt = ["17V.CSV","18V.CSV",
         "19V.CSV","20V.CSV","21V.CSV",
         "22V.CSV","23V.CSV","24V.CSV","25V.CSV",
         "26V.CSV","27V.CSV","28V.CSV","29V.CSV",
         "30V.CSV","31V.CSV","32V.CSV","33V.CSV",
         "34V.CSV","35V.CSV","36V.CSV","37V.CSV",
         "38V.CSV","39V.CSV","40V.CSV","41V.CSV",
         "42V.CSV","43V.CSV","44V.CSV","45V.CSV",
         "46V.CSV","47V.CSV","48V.CSV",
         "49V.CSV","50V.CSV"]
keys_volt = ["17V","18V",
         "19V","20V","21V",
         "22V","23V","24V","25V",
         "26V","27V","28V","29V",
         "30V","31V","32V","33V",
         "34V","35V","36V","37V",
         "38V","39V","40V","41V",
         "42V","43V","44V","45V",
         "46V","47V","48V",
         "49V","50V"]
voltage_volt = np.array([17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,
                    40,41,42,43,44,45,46,47,48,49,50])

counter=17
data={}
for i in files_volt:
    data["{0}V".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=1


# max_volt=np.array([])
t_0_volt=np.array([])
unc_t_0_volt=np.array([])
fwhm_volt=np.array([])
unc_fwhm_volt=np.array([])
n_0_volt=np.array([])
unc_n_0_volt=np.array([])

count=0 
counter=17
for key in keys_volt:

    if 17<=counter<=22:

        time = data[key][500:2400,3]
        volt = data[key][500:2400,4] 
        peaks, properties = find_peaks(volt, distance=1500, width=100, height = -8)

        volt_left = volt[:int(properties["left_ips"])-250]
        volt_right = volt[int(properties["right_ips"])+250:]
        volt_middle = np.full(int(properties["right_ips"])+250-(int(properties["left_ips"])-250), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time
    
    if 22<counter<35:
        time = data[key][400:2000,3]
        volt = data[key][400:2000,4]

        peaks, properties = find_peaks(volt, distance=1000, width=50, height = -8)

    
        volt_left = volt[:int(properties["left_ips"])-200]
        volt_right = volt[int(properties["right_ips"])+200:]
        volt_middle = np.full(int(properties["right_ips"])+200-(int(properties["left_ips"])-200), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time  
        
    if counter>34:
        time = data[key][400:1600,3]
        volt = data[key][400:1600,4]   

        peaks, properties = find_peaks(volt, distance=1000, width=50, height = -8)
        volt_left = volt[:int(properties["left_ips"])-200]
        volt_right = volt[int(properties["right_ips"])+200:]
        volt_middle = np.full(int(properties["right_ips"])+200-(int(properties["left_ips"])-200), np.average((volt_right[0],volt_left[-1])))
        volt_concat = np.concatenate((volt_left,volt_middle,volt_right))
        time_concat = time          
    
        
        
    popt_exp, pcov_exp = curve_fit(exponential, time_concat, volt_concat, p0=(10, -1, -1.5),check_finite=False, maxfev=10000000)
    background = exponential(time, popt_exp[0],popt_exp[1],popt_exp[2])
    volt_subt = volt - (background-background[-1])     
    
    popt_gau, pcov_gau = curve_fit(gaussian, time, volt_subt,  p0=(7, 2e-5, 1e-5, 8))
    unc_popt_gau = np.sqrt(np.diag(pcov_gau))

    t_0_volt = np.append(t_0_volt, popt_gau[1])
    n_0_volt = np.append(n_0_volt, popt_gau[0])        
    fwhm_volt = np.append(fwhm_volt, popt_gau[2]*2*np.sqrt(2*np.log(2)))
    
    unc_t_0_volt = np.append(unc_t_0_volt, unc_popt_gau[1])
    unc_n_0_volt = np.append(unc_n_0_volt, unc_popt_gau[0])
    unc_fwhm_volt = np.append(unc_fwhm_volt, unc_popt_gau[2]*2*np.sqrt(2*np.log(2)))
    
    
    # plt.scatter(time, volt_subt-popt_gau[3], marker='.', color="r")
    # plt.plot(time, gaussian(time, popt_gau[0], popt_gau[1], popt_gau[2], popt_gau[3])-popt_gau[3])
    # plt.xlabel("time [s]")
    # plt.ylabel("voltage [V]")
    # plt.title("Signal for d=0.3cm and V$_s$={0}V".format(counter))
    # plt.savefig("gaussian{0}".format(counter),dpi=300)
    # plt.show()

    count+=1
    counter+=1


files_d = ["3.0mm.CSV","3.1mm.CSV",
         "3.2mm.CSV","3.3mm.CSV","3.4mm.CSV",
         "3.5mm.CSV","3.6mm.CSV","3.7mm.CSV","3.8mm.CSV",
         "3.9mm.CSV","4.0mm.CSV","4.1mm.CSV","4.2mm.CSV",
         "4.3mm.CSV","4.4mm.CSV","4.5mm.CSV","4.6mm.CSV",
         "4.7mm.CSV","4.8mm.CSV","4.9mm.CSV","5.0mm.CSV",
         "5.1mm.CSV","5.2mm.CSV","5.3mm.CSV","5.4mm.CSV",
         "5.5mm.CSV"]
keys_d = ["30","31","32","33","34",
         "35","36","37","38",
         "39","40","41","42",
         "43","44","45","46",
         "47","48","49","50",
         "51","52","53","54",
         "55"]
d_d = np.array([3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.5])
unc_d_d=0.05

counter=30
data={}
for i in files_d:
    data["{0}".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=1



t_0_d=np.array([])
unc_t_0_d=np.array([])
fwhm_d=np.array([])
unc_fwhm_d=np.array([])
n_0_d=np.array([])
unc_n_0_d=np.array([])

count=0 
counter=30
for key in keys_d:

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

    t_0_d = np.append(t_0_d, popt_gau[1])
    n_0_d = np.append(n_0_d, popt_gau[0])        
    fwhm_d = np.append(fwhm_d, popt_gau[2]*2*np.sqrt(2*np.log(2)))
    
    unc_t_0_d = np.append(unc_t_0_d, unc_popt_gau[1])
    unc_n_0_d = np.append(unc_n_0_d, unc_popt_gau[0])
    unc_fwhm_d = np.append(unc_fwhm_d, unc_popt_gau[2]*2*np.sqrt(2*np.log(2)))
    
    # plt.scatter(time, volt_subt-popt_gau[3], marker='.', color="r")
    # plt.plot(time, gaussian(time, popt_gau[0], popt_gau[1], popt_gau[2], popt_gau[3])-popt_gau[3])
    # plt.xlabel("time [s]")
    # plt.ylabel("voltage [V]")
    # plt.title("Signal for V$_s$=30V and d={:.1f}mm".format(counter*0.1))
    # plt.savefig("gaussian{0}".format(counter),dpi=300)
    # plt.show()

    count+=1
    counter+=1
 
l=2.0 #in cm
unc_l=0.1 #in cm
voltage_d = 30.0
unc_voltage_d = 0.1

d_d*=0.1 #in cm
unc_d_d*=0.1 # in cm

log_d = np.log(np.sqrt(t_0_d)*abs(n_0_d))
    

d_volt=0.3 #in cm
unc_d_volt=0.005 #in cm

unc_voltage_volt=0.1 # in V

log_volt = np.log(np.sqrt(t_0_volt)*abs(n_0_volt))





popt_mo_volt, pcov_mo_volt = curve_fit(proportional, voltage_volt/l, d_volt/t_0_volt)
unc_popt_mo_volt = np.sqrt(np.diag(pcov_mo_volt))
mobility_volt = popt_mo_volt[0]
unc_mobility_volt = unc_popt_mo_volt[0]
print("Mobility = ", mobility_volt, "+/-", unc_mobility_volt, "cm^2 s^-1 V-1, rel = ", unc_mobility_volt/mobility_volt)


popt_mo_d, pcov_mo_d = curve_fit(proportional, voltage_d/d_d, l/t_0_d)
unc_popt_mo_d = np.sqrt(np.diag(pcov_mo_d))
mobility_d = popt_mo_d[0]
unc_mobility_d = unc_popt_mo_d[0]
print("Mobility = ", mobility_d, "+/-", unc_mobility_d, "cm^2 s^-1 V-1, rel = ", unc_mobility_d/mobility_d)



popt_dif_volt, pcov_dif_volt = curve_fit(proportional, t_0_volt**3, (fwhm_volt*d_volt)**2)
unc_popt_dif_volt = np.sqrt(np.diag(pcov_dif_volt))
diffusivity_volt = popt_dif_volt[0]/(16*np.log(2))
unc_diffusivity_volt = diffusivity_volt*np.sqrt((unc_popt_dif_volt[0]/popt_dif_volt[0])**2+(2*unc_d_volt/d_volt)**2)
print("Diffusivity = ", diffusivity_volt, "+/-", unc_diffusivity_volt, "cm^2 s^-1, rel = ", unc_diffusivity_volt/diffusivity_volt)


popt_dif_d, pcov_dif_d = curve_fit(proportional, t_0_d**3, (d_d*fwhm_d)**2)
unc_popt_dif_d = np.sqrt(np.diag(pcov_dif_d))
diffusivity_d = popt_dif_d[0]/(16*np.log(2))
unc_diffusivity_d = diffusivity_d*(unc_popt_dif_d[0]/popt_dif_d[0])
print("Diffusivity = ", diffusivity_d, "+/-", unc_diffusivity_d, "cm^2 s^-1, rel = ", unc_diffusivity_d/diffusivity_d)




x_data = np.linspace(0, 3e-14, 1000)

# y_unc_volt = 2*fwhm_volt*unc_fwhm_volt
plt.plot(x_data, proportional(x_data, popt_dif_volt[0]))
plt.scatter(t_0_volt**3, (fwhm_volt*d_volt)**2, marker=".", label="varying V")
#plt.errorbar(t_0**3, fwhm**2, yerr=y_unc, ls="none")

plt.plot(x_data, proportional(x_data, popt_dif_d[0]))
plt.scatter(t_0_d**3, (fwhm_d*d_d)**2, marker=".", label="varying d")
plt.legend()

plt.xlabel(r"$t_0^3$ [s$^3$]")
plt.ylabel(r"$(t_p d)^2$ [s$^3$ cm$^2$]")
plt.title("Diffusion Constant")
plt.savefig("diffusion_fpb",dpi=300)
plt.show()



# non linear lifetime
popt_nlin_volt, pcov_nlin_volt = curve_fit(non_linear, t_0_volt, abs(n_0_volt), p0=(0.02,-100000))
unc_popt_nlin_volt = np.sqrt(np.diag(pcov_nlin_volt))
lifetime_volt = -1/popt_nlin_volt[1]
unc_lifetime_volt = lifetime_volt*abs(unc_popt_nlin_volt[1]/popt_nlin_volt[1])
print("Lifetime (non linear) = ", lifetime_volt, "+/-", unc_lifetime_volt, "s, rel = ", unc_lifetime_volt/lifetime_volt)

# x_data = np.linspace(1.1e-5,2.6e-5,1000)
# plt.scatter(t_0,abs(n_0),marker='.')
# plt.plot(x_data, non_linear(x_data, popt_nlin[0], popt_nlin[1]))
# plt.errorbar(t_0, abs(n_0), yerr=unc_n_0,ls="none")
# plt.xlabel(r"$t_0$ [s]")
# plt.ylabel(r"$n_0$")
# plt.title("Lifetime for d=0.3cm varying $V_s$")
# plt.savefig("lifetime_non_linear", dpi=300)
# plt.show()

# mobility correction analysis
mobility_array_volt = (d_volt*l) / (t_0_volt*voltage_volt)
aveg_mob = np.average(mobility_array_volt)
x = (2*constants.k*300*l) / (constants.e*voltage_volt*d_volt) * (t_0_volt / lifetime_volt + 0.5)
mobility_correction = mobility_array_volt * (np.sqrt(1 + x**2) - x)
plt.scatter(voltage_volt, mobility_array_volt, marker=".", label="measured values")
plt.scatter(voltage_volt, mobility_correction, marker=".", label="corrected values")

# mobility_array_d = (d_d*l) / (t_0_d*voltage_d)
# #mobility_array_unc_d = mobility_array_d*np.sqrt((unc_d_d/d_d)**2+(unc_l/l)**2+(unc_t_0_d/t_0_d)**2+(unc_voltage_d/voltage_d)**2)
# plt.scatter(d_d, mobility_array_d, marker=".")
# #plt.errorbar(d_d,mobility_array,yerr=mobility_array_unc, ls='none')


plt.legend()
plt.ylim([0,1550])
plt.xlabel(r"$V_s$ [V]")
plt.ylabel(r"$\mu$ [cm$^2$ s$^{-1}$ V$^{-1}$]")
plt.title(r"Mobility $")
plt.savefig("mobility_horizontal_fpb", dpi=300)
plt.show()




# mobility_corrected = np.average(mobility_correction)
# mobility_corrected_unc = np.std(mobility_correction)
# print("Mobility final = ", mobility_corrected, "+/-", mobility_corrected_unc, "cm^2 s^-1 V-1, rel = ", mobility_corrected_unc/mobility_corrected)

# # einstein ratio
# einstein_ratio = diffusivity / mobility_corrected
# print("Einstein relation = ", einstein_ratio, "V")

