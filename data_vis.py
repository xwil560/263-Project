from matplotlib import pyplot as plt
import numpy as np
import os

def pilot_plots():

    # Creating plots for pilot study data
    plt.rcParams["figure.figsize"] = (8,6)
    os.chdir("data")
    f, axs = plt.subplots(3, sharex=True)

    # Importing steam injection data
    steam_time = np.genfromtxt('tr_steam.txt',skip_header=True,delimiter=",",usecols=0)
    steam = np.genfromtxt('tr_steam.txt',skip_header=True,delimiter=",",usecols=1)
    axs[0].plot(steam_time,steam,'k-',label="Steam Rate (tonnes/day)")

    # Importing oil and water extraction data onto same axis
    oil_time = np.genfromtxt('tr_oil.txt',skip_header=True,delimiter=",",usecols=0)
    oil = np.genfromtxt('tr_oil.txt',skip_header=True,delimiter=",",usecols=1)
    axs[1].plot(oil_time,oil,'b-',label="Oil Rate ($m^{3}$/day)")
    axs4 = axs[1].twinx()
    water_time = np.genfromtxt('tr_water.txt',skip_header=True,delimiter=",",usecols=0)
    water = np.genfromtxt('tr_water.txt',skip_header=True,delimiter=",",usecols=1)
    axs4.plot(water_time,water,'m-',label="Water Rate ($m^{3}$/day)")

    # Importing temperature and pressure data onto same axis
    temp_time = np.genfromtxt('tr_T.txt',skip_header=True,delimiter=",",usecols=0)
    temp = np.genfromtxt('tr_T.txt',skip_header=True,delimiter=",",usecols=1)
    axs[2].plot(temp_time,temp,'r-',label="Temperature ($^{0}C$)")
    axs3 = axs[2].twinx()
    pressure_time = np.genfromtxt('tr_p.txt',skip_header=True,delimiter=",",usecols=0)
    pressure = np.genfromtxt('tr_p.txt',skip_header=True,delimiter=",",usecols=1)
    axs3.plot(pressure_time,pressure,'g-',label="Pressure (kPa)")

    # Labelling overall plot title
    plt.suptitle("Pilot Study Time Series")

    # Labelling axis and labels, as well as general formatting of the legends and axis to create easier visualisation
    axs[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[0].set_ylabel("Steam Rate (tonnes/day)")
    axs[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[1].set_ylabel("Oil Rate ($m^{3}$/day)")
    axs[1].tick_params(axis='y',colors='blue')
    axs[1].yaxis.label.set_color('blue')
    axs[2].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[2].set_xlabel("Time (days)")
    axs[2].set_ylabel("Temperature ($^{0}C$)")
    axs[2].tick_params(axis='y',colors='red')
    axs[2].yaxis.label.set_color('red')
    axs3.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',ncol=2, borderaxespad=0.)
    axs3.set_ylabel("Pressure (kPa)")
    axs3.tick_params(axis='y',colors='green')
    axs3.yaxis.label.set_color('green')
    axs4.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',ncol=2, borderaxespad=0.)
    axs4.set_ylabel("Water Rate ($m^{3}$/day)")
    axs4.tick_params(axis='y',colors='purple')
    axs4.yaxis.label.set_color('purple')

    # Changing subdirectory, ensuring tight layout for plot, and both saving and showing the plot
    os.chdir("../plots")
    plt.tight_layout()
    plt.savefig('pilot_ts.png', dpi=300)
    plt.show()
    os.chdir("..")

def exp2_plots():
    
    # Creating plots for experiment 2 data
    plt.rcParams["figure.figsize"] = (8,6)
    os.chdir("pilot_project_data_2012")
    f, axs = plt.subplots(3, sharex=True)

    # Importing steam injection data
    steam_time = np.genfromtxt('exp2_steam.csv',skip_header=True,delimiter=",",usecols=0)
    steam = np.genfromtxt('exp2_steam.csv',skip_header=True,delimiter=",",usecols=1)
    axs[0].plot(steam_time,steam,'k-',label="Steam Rate (tonnes/day)")

    # Importing oil and water extraction data onto same axis
    oil_time = np.genfromtxt('exp2_oil.csv',skip_header=True,delimiter=",",usecols=0)
    oil = np.genfromtxt('exp2_oil.csv',skip_header=True,delimiter=",",usecols=1)
    axs[1].plot(oil_time,oil,'b-',label="Oil Rate ($m^{3}$/day)")
    axs4 = axs[1].twinx()
    water_time = np.genfromtxt('exp2_water.csv',skip_header=True,delimiter=",",usecols=0)
    water = np.genfromtxt('exp2_water.csv',skip_header=True,delimiter=",",usecols=1)
    axs4.plot(water_time,water,'m-',label="Water Rate ($m^{3}$/day)")

    # Importing temperature and pressure data onto same axis
    temp_time = np.genfromtxt('exp2_temp.csv',skip_header=True,delimiter=",",usecols=0)
    temp = np.genfromtxt('exp2_temp.csv',skip_header=True,delimiter=",",usecols=1)
    axs[2].plot(temp_time,temp,'r-',label="Temperature ($^{0}C$)")
    axs3 = axs[2].twinx()
    pressure_time = np.genfromtxt('exp2_pressure.csv',skip_header=True,delimiter=",",usecols=0)
    pressure = np.genfromtxt('exp2_pressure.csv',skip_header=True,delimiter=",",usecols=1)
    axs3.plot(pressure_time,pressure,'g-',label="Pressure (kPa)")

    # Labelling overall plot title
    plt.suptitle("Experiment 2 Time Series")

    # Labelling axis and labels, as well as general formatting of the legends and axis to create easier visualisation
    axs[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[0].set_ylabel("Steam Rate (tonnes/day)")
    axs[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[1].set_ylabel("Oil Rate ($m^{3}$/day)")
    axs[1].tick_params(axis='y',colors='blue')
    axs[1].yaxis.label.set_color('blue')
    axs[2].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[2].set_xlabel("Time (days)")
    axs[2].set_ylabel("Temperature ($^{0}C$)")
    axs[2].tick_params(axis='y',colors='red')
    axs[2].yaxis.label.set_color('red')
    axs3.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',ncol=2, borderaxespad=0.)
    axs3.set_ylabel("Pressure (kPa)")
    axs3.tick_params(axis='y',colors='green')
    axs3.yaxis.label.set_color('green')
    axs4.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',ncol=2, borderaxespad=0.)
    axs4.set_ylabel("Water Rate ($m^{3}$/day)")
    axs4.tick_params(axis='y',colors='purple')
    axs4.yaxis.label.set_color('purple')

    # Changing subdirectory, ensuring tight layout for plot, and both saving and showing the plot
    os.chdir("../plots")
    plt.tight_layout()
    plt.savefig('exp2_ts.png', dpi=300)
    plt.show()
    os.chdir("..")

def exp3_plots():
    
    # Creating plots for experiment 3 data
    plt.rcParams["figure.figsize"] = (8,6)
    os.chdir("pilot_project_data_2012")
    f, axs = plt.subplots(3, sharex=True)

    # Importing steam injection data
    steam_time = np.genfromtxt('exp3_steam.csv',skip_header=6,delimiter=",",usecols=0)
    steam = np.genfromtxt('exp3_steam.csv',skip_header=6,delimiter=",",usecols=1)
    axs[0].plot(steam_time,steam,'k-',label="Steam Rate (tonnes/day)")

    # Importing oil and water extraction data onto same axis
    oil_time = np.genfromtxt('exp3_oil.csv',skip_header=6,delimiter=",",usecols=0)
    oil = np.genfromtxt('exp3_oil.csv',skip_header=6,delimiter=",",usecols=1)
    axs[1].plot(oil_time,oil,'b-',label="Oil Rate ($m^{3}$/day)")
    axs4 = axs[1].twinx()
    water_time = np.genfromtxt('exp3_water.csv',skip_header=6,delimiter=",",usecols=0)
    water = np.genfromtxt('exp3_water.csv',skip_header=6,delimiter=",",usecols=1)
    axs4.plot(water_time,water,'m-',label="Water Rate ($m^{3}$/day)")

    # Importing temperature and pressure data onto same axis
    temp_time = np.genfromtxt('exp3_temp.csv',skip_header=6,delimiter=",",usecols=0)
    temp = np.genfromtxt('exp3_temp.csv',skip_header=6,delimiter=",",usecols=1)
    axs[2].plot(temp_time,temp,'r-',label="Temperature ($^{0}C$)")
    axs3 = axs[2].twinx()
    pressure_time = np.genfromtxt('exp3_pressure.csv',skip_header=6,delimiter=",",usecols=0)
    pressure = np.genfromtxt('exp3_pressure.csv',skip_header=6,delimiter=",",usecols=1)
    axs3.plot(pressure_time,pressure,'g-',label="Pressure (kPa)")

    # Labelling overall plot title
    plt.suptitle("Experiment 3 Time Series")

    # Labelling axis and labels, as well as general formatting of the legends and axis to create easier visualisation
    axs[0].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[0].set_ylabel("Steam Rate (tonnes/day)")
    axs[1].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[1].set_ylabel("Oil Rate ($m^{3}$/day)")
    axs[1].tick_params(axis='y',colors='blue')
    axs[1].yaxis.label.set_color('blue')
    axs[2].legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower left',ncol=2, borderaxespad=0.)
    axs[2].set_xlabel("Time (days)")
    axs[2].set_ylabel("Temperature ($^{0}C$)")
    axs[2].tick_params(axis='y',colors='red')
    axs[2].yaxis.label.set_color('red')
    axs3.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',ncol=2, borderaxespad=0.)
    axs3.set_ylabel("Pressure (kPa)")
    axs3.tick_params(axis='y',colors='green')
    axs3.yaxis.label.set_color('green')
    axs4.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc='lower right',ncol=2, borderaxespad=0.)
    axs4.set_ylabel("Water Rate ($m^{3}$/day)")
    axs4.tick_params(axis='y',colors='purple')
    axs4.yaxis.label.set_color('purple')

    # Changing subdirectory, ensuring tight layout for plot, and both saving and showing the plot
    os.chdir("../plots")
    plt.tight_layout()
    plt.savefig('exp3_ts.png', dpi=300)
    plt.show()
    os.chdir("..")