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
#"12V.CSV","13V.CSV","14V.CSV","15V.CSV","16V.CSV","17V.CSV","18V.CSV",
#"12V","13V","14V","15V","16V","17V","18V",
#12,13,14,15,16,17,18,
files = ["15V.CSV","16V.CSV","17V.CSV","18V.CSV",
         "19V.CSV","20V.CSV","21V.CSV",
         "22V.CSV","23V.CSV","24V.CSV","25V.CSV",
         "26V.CSV","27V.CSV","28V.CSV","29V.CSV"]
keys = ["15V","16V","17V","18V",
         "19V","20V","21V",
         "22V","23V","24V","25V",
         "26V","27V","28V","29V"]
voltage = np.array([15,16,17,18,19,20,21,22,23,24,25,26,27,28,29])


counter=15
data={}
for i in files:
    data["{0}V".format(counter)]=np.genfromtxt(i, dtype=float, delimiter=",")
    counter+=1

max_volt=np.array([])
t_0=np.array([])
unc_t_0=np.array([])

count=0
for key in keys:
    
    if count<2:
        time=data[key][1200:2000,3]
        volt=data[key][1200:2000,4]
        print(max(volt))
    else:
        time=data[key][500:2200,3]
        volt=data[key][500:2200,4]
    
    max_volt=np.append(max_volt,max(volt))
    indexes=np.where(volt==max_volt[count])
    count+=1
    
    if len(indexes[0])==1:
        t_0=np.append(t_0,time[indexes[0]])
        unc_t_0= np.append(unc_t_0,(time[1]-time[0]))
    
    else:
        time_range = np.array([])
        last_index=indexes[0][-1]
        for i in range(indexes[0][0],last_index,1):
            time_range=np.append(time_range, time[i])
        t_0 = np.append(t_0,np.mean(time_range)) 
        unc_t_0= np.append(unc_t_0,time_range[-1]-time_range[0])
    
 
rel_unc=unc_t_0/t_0


x_axis = voltage/(20*0.1)
y_axis = 3.0*0.1/t_0
y_unc = y_axis*np.sqrt((unc_t_0/t_0)**2+(0.1/3.0)**2)
fit = curve_fit(linear, x_axis, y_axis,p0=(500,))

x_plot = np.linspace(0,15,100)
mobility=fit[0][0]
unc_mobility = np.sqrt(np.diag(fit[1]))
print("mobility= ", mobility, "+/-", unc_mobility)

chi2=chi_squared(y_axis, linear(fit[0][0],x_axis), y_unc)
red_chi2=chi2/(len(x_axis)-1)

print(red_chi2)

plt.scatter(x_axis, y_axis,marker=".")
plt.errorbar(x_axis,y_axis, yerr=y_unc,ls='none')
plt.plot(x_plot, fit[0][0]*x_plot)
plt.xlabel("V_s/l")
plt.ylabel("d/t")
plt.savefig("mobility_day3",dpi=300)
plt.show()

