# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import norm


def linear(m,x):
    return m*x

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

files = ["17V.CSV","19V.CSV","21V.CSV","23V.CSV","25V.CSV","27V.CSV","29V.CSV"]
keys = ["17V","19V","21V","23V","25V","27V","29V"]
voltage = np.array([17,19,21,23,25,27,29])



counter=17
data={}
time = np.array([])
voltage = np.array([])
for i in files:
    data["{0}V".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=2


max_volt=np.array([])
t_0=np.array([])
unc_t_0=np.array([])
fwhm=np.array([])
unc_fwhm=np.array([])
chi2 = np.array([])
count=0 
counter=17   
for key in keys:

    time=data[key][800:2450,3]
    volt=data[key][800:2450,4]*5e-2
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
    unc_popt_gau = np.sqrt(np.diag(pcov_gau))
    chi2 = np.append(chi2, chi_squared(volt, gaussian(time, popt_gau[0],popt_gau[1],popt_gau[2],popt_gau[3]), unc_volt)/len(volt))

    fwhm = np.append(fwhm, 2*np.sqrt(2*np.log(2))*popt_gau[2])
    #unc_fwhm = np.append(unc_fwhm, fwhm[count]*unc_popt_gau[2]/popt_gau[2])
    #unc_fwhm = np.append(unc_fwhm, 2*np.sqrt(2*np.log(2))*unc_popt_gau[2])
    #print(fwhm[count], "+/-", unc_fwhm[count])
    
    plt.plot(time, gaussian(time, popt_gau[0],popt_gau[1],popt_gau[2],popt_gau[3]), marker=".")
    plt.scatter(time, volt, marker='.', color="r")
    plt.xlabel("time")
    plt.ylabel("voltage")
    plt.savefig("gaussian{0}".format(counter),dpi=300)
    plt.show()
    
    count+=1
    counter+=2


popt, pcov = curve_fit(linear,t_0**3, fwhm**2)
unc_popt = np.sqrt(np.diag(pcov))

x_data = np.linspace(0, 1.7e-14, 1000)
plt.plot(x_data, linear(popt[0], x_data))

plt.scatter(t_0**3, fwhm**2, marker=".")
plt.errorbar(t_0**3, fwhm**2, yerr=2*fwhm**2*unc_popt_gau[2]/popt_gau[2], ls="none")
plt.xlabel("t_0^3")
plt.ylabel("t_p^2")
plt.savefig("diffusion constant",dpi=300)
plt.show()


d=0.3 #in cm
unc_d=0.005 #in cm

dif = popt[0]*0.3**2/(16*np.log(2))
unc_dif=dif*np.sqrt((unc_popt[0]/popt[0])**2+(2*unc_d/d)**2)
print("diffusion constant=", dif, "+/-", unc_dif, ", ", unc_dif/dif, "%")

