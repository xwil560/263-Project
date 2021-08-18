# imports
import numpy as np
from matplotlib import pyplot as plt
#from .py import *


def temperature_benchmark():
    ''' Check if the temperature fit the model

        Parameter:
        -------
        None

        Return
        --------
        None
    '''
    # solve the ode numerically
    #ts,xsn=solve_ode(ode_model, t0=0, t1=10, dt=1/10,x0=0, pars=[-1,1,1,0])
    # solve the ode with the given parameter analytically
    # xsa=np.exp(-ts)-1
    # subplot the 2 graphs
    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(ts, xsn, '-x')
    ax1.plot(ts, xsa, '-')
    ax2.plot(ts, np.absolute(xsn-xsa)/np.absolute(xsa), '-')
    dtc = np.linspace(1, 3, 21)
    lastval = 0*dtc
    for i in range(21):
        tsc, xsnc = solve_ode(ode_model, t0=0, t1=10,
                              dt=1/dtc[i], x0=0, pars=[-1, 1, 1, 0])
        lastval[i] = xsnc[-1]
    ax3.plot(dtc, lastval, 'x')
    plt.show()


def pressure_benchmark():
    ''' Check if the pressure fit the model

        Parameter:
        ---------
        None

        Return:
        -------
        None
    '''
    # solve the ode numerically
    #ts,xsn=solve_ode(ode_model, t0=0, t1=10, dt=1/10,x0=0, pars=[-1,1,1,0])
    # solve the ode with the given parameter analytically
    # xsa=np.exp(-ts)-1
    # subplot the 2 graphs
    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(ts, xsn, '-x')
    ax1.plot(ts, xsa, '-')
    ax2.plot(ts, np.absolute(xsn-xsa)/np.absolute(xsa), '-')
    dtc = np.linspace(1, 3, 21)
    lastval = 0*dtc
    for i in range(21):
        tsc, xsnc = solve_ode(ode_model, t0=0, t1=10,
                              dt=1/dtc[i], x0=0, pars=[-1, 1, 1, 0])
        lastval[i] = xsnc[-1]
    ax3.plot(dtc, lastval, 'x')
    plt.show()
