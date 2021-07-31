from matplotlib import pyplot as plt

import numpy as np
import os
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
os.chdir("data")
f, ax1 = plt.subplots(1, 1)
time1=np.genfromtxt('tr_oil.txt',skip_header=True,delimiter=",",usecols=0)
vol1=np.genfromtxt('tr_oil.txt',skip_header=True,delimiter=",",usecols=1)
ax1.plot(time1,vol1,'-',label="Oil extraction rate ($m^{3}$/days)")

os.chdir("..")


os.chdir("data")

time2=np.genfromtxt('tr_water.txt',skip_header=True,delimiter=",",usecols=0)
vol2=np.genfromtxt('tr_water.txt',skip_header=True,delimiter=",",usecols=1)
ax1.plot(time2,vol2,'-',label="Water extraction rate ($m^{3}$/days)")

os.chdir("..")


os.chdir("data")

time3=np.genfromtxt('tr_steam.txt',skip_header=True,delimiter=",",usecols=0)
mass3=np.genfromtxt('tr_steam.txt',skip_header=True,delimiter=",",usecols=1)
ax1.plot(time3,mass3,'-',label="Steam injection rate (tonnes/day)")

os.chdir("..")


os.chdir("data")

time4=np.genfromtxt('tr_T.txt',skip_header=True,delimiter=",",usecols=0)
temp=np.genfromtxt('tr_T.txt',skip_header=True,delimiter=",",usecols=1)
ax1.plot(time4,temp,'-',label="Temperature at 350m ($^{0}C$)")

os.chdir("..")


os.chdir("data")

time=np.genfromtxt('tr_p.txt',skip_header=True,delimiter=",",usecols=0)
pressure=np.genfromtxt('tr_p.txt',skip_header=True,delimiter=",",usecols=1)
ax1.plot(time,pressure,'-',label="Pressure at 350m (kPa)")
plt.title("Data visulisation")
plt.xlabel("Time (days)")
plt.legend(bbox_to_anchor=(1.05,1),fontsize='xx-small')
os.chdir("..")
plt.savefig('data.png')
