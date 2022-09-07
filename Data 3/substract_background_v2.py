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
            
    plt.scatter(time, volt, marker='.', color="b")
    plt.plot(time, exponential(time, popt_exp[0],popt_exp[1],popt_exp[2]), marker=".")
    plt.xlabel("time")
    plt.ylabel("voltage")
    plt.title("{0}V".format(counter))
    plt.savefig("gaussian{0}_exp_background".format(counter),dpi=300)
    plt.show()
    
    background = exponential(time, popt_exp[0],popt_exp[1],popt_exp[2])
    volt_subt = volt - (background-background[-1]) 
    plt.scatter(time, volt_subt, marker='.', color="r")
    plt.xlabel("time")
    plt.ylabel("voltage")
    plt.title("{0}V".format(counter))
    plt.show()


    count+=1
    counter+=1