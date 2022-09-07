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

files = ["15V.CSV","20V.CSV","25V.CSV","30V.CSV","35V.CSV","40V.CSV","45V.CSV"]
keys = ["15V","20V","25V","30V","35V","40V","45V"]
voltage = np.array([15,20,25,30,35,40,45])

counter=15
data={}
for i in files:
    data["{0}V".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=5

max_volt=np.array([])
t_0=np.array([])
unc_t_0=np.array([])

count=0
for key in keys:

    #indexes=0
    times=0
    volt=0
    last_index=0

    
    time=data[key][1000:2200,3]
    volt=data[key][1000:2200,4]
    
   

    max_volt=np.append(max_volt,max(volt))
    indexes=np.where(volt==max_volt[count])
    count+=1
    
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
    
    
    
rel_unc=unc_t_0/t_0


x_axis = voltage/(20*0.1)
y_axis = 1.79*0.1/t_0
y_unc = y_axis*np.sqrt((unc_t_0/t_0)**2+(0.05/1.79)**2)
fit = curve_fit(linear, x_axis, y_axis,p0=(1000,))

x_plot = np.linspace(0,23,100)
print("mobility= ", fit[0][0])

chi2=chi_squared(y_axis, linear(fit[0][0],x_axis), y_unc)
red_chi2=chi2/len(x_axis)

print(red_chi2)

plt.scatter(x_axis, y_axis,marker=".")
plt.errorbar(x_axis,y_axis, yerr=y_unc,ls='none')
plt.plot(x_plot, fit[0][0]*x_plot)
plt.xlabel("V_s/l")
plt.ylabel("d/t")
plt.savefig("mobility_v1",dpi=300)
plt.show()



