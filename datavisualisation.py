from matplotlib import pyplot as plt
import numpy as np
import os
oil=True
water=True
steam=True
temp=True
pressure=True
#plotting oil production rate 
if oil:
    os.chdir("data")
    f, ax1 = plt.subplots(1, 1)
    time=np.genfromtxt('tr_oil.txt',skip_header=True,delimiter=",",usecols=0)
    vol=np.genfromtxt('tr_oil.txt',skip_header=True,delimiter=",",usecols=1)
    ax1.plot(time,vol,'-')
    ax1.set_title('Oil production rate')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Production rate (m^3/day)')
    os.chdir("..")
    plt.savefig('oil.png')
if water:
    os.chdir("data")
    f, ax1 = plt.subplots(1, 1)
    time=np.genfromtxt('tr_water.txt',skip_header=True,delimiter=",",usecols=0)
    vol=np.genfromtxt('tr_water.txt',skip_header=True,delimiter=",",usecols=1)
    ax1.plot(time,vol,'-')
    ax1.set_title('Water production rate')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Production rate (m^3/day)')
    os.chdir("..")
    plt.savefig('water.png')
if steam:
    os.chdir("data")
    f, ax1 = plt.subplots(1, 1)
    time=np.genfromtxt('tr_steam.txt',skip_header=True,delimiter=",",usecols=0)
    mass=np.genfromtxt('tr_steam.txt',skip_header=True,delimiter=",",usecols=1)
    ax1.plot(time,mass,'-')
    ax1.set_title('Steam injection rate')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Injection rate (tonnes/day)')
    os.chdir("..")
    plt.savefig('steam.png')
if temp:
    os.chdir("data")
    f, ax1 = plt.subplots(1, 1)
    time=np.genfromtxt('tr_T.txt',skip_header=True,delimiter=",",usecols=0)
    temp=np.genfromtxt('tr_T.txt',skip_header=True,delimiter=",",usecols=1)
    ax1.plot(time,temp,'-')
    ax1.set_title('Temperature at 350m into the well')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Temperature (degree C)')
    os.chdir("..")
    plt.savefig('temp.png')
if pressure:
    os.chdir("data")
    f, ax1 = plt.subplots(1, 1)
    time=np.genfromtxt('tr_p.txt',skip_header=True,delimiter=",",usecols=0)
    temp=np.genfromtxt('tr_p.txt',skip_header=True,delimiter=",",usecols=1)
    ax1.plot(time,temp,'-')
    ax1.set_title('Pressure at 350m into the well')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Pressure (kPa)')
    os.chdir("..")
    plt.savefig('pressure.png')
