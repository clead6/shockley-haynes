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
counter=16
for key in keys:
    
    if count<2:
        time=data[key][1000:2000,3]
        volt=data[key][1000:2000,4]

    else:
        time=data[key][500:2000,3]
        volt=data[key][500:2000,4]
    

    #unc_volt = 0.0005
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


    popt_gau, pcov_gau = curve_fit(gaussian, time, volt, p0=(7, 2e-5, 1e-5, 8))
    unc_popt_gau = np.sqrt(np.diag(pcov_gau))
    chi2 = np.append(chi2, chi_squared(volt, gaussian(time, popt_gau[0],popt_gau[1],popt_gau[2],popt_gau[3]), unc_volt)/len(volt))
    unc_t_0= np.append(unc_t_0, abs(popt_gau[2]) / 2 )
    
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
    counter+=1


x_axis = voltage/(20*0.1)
y_axis = 3.0*0.1/t_0
print(unc_t_0/t_0)
y_unc = y_axis*np.sqrt((unc_t_0/t_0)**2+(0.1/3.0)**2)
fit = curve_fit(linear, x_axis, y_axis,p0=(500,))

x_plot = np.linspace(0,15,100)
mobility=fit[0][0]
unc_mobility = np.sqrt(np.diag(fit[1]))
print("mobility= ", mobility, "+/-", unc_mobility[0], ", ", unc_mobility[0]/mobility *100, "%")

chi2=chi_squared(y_axis, linear(fit[0][0],x_axis), y_unc)
red_chi2=chi2/(len(x_axis)-1)

print(red_chi2)

plt.scatter(x_axis, y_axis,marker=".")
plt.errorbar(x_axis,y_axis, yerr=y_unc,ls='none')
plt.plot(x_plot, fit[0][0]*x_plot)
plt.xlabel("V_s/l")
# plt.ylabel("d/t")
plt.savefig("mobility_day3_different_unc",dpi=300)
plt.show()

