# imports
import numpy as np
from matplotlib import pyplot as plt
from LPM import *


def pressure_benchmark():
    ''' Check if the pressure ode is implemented correctly by comparing it to the analytical solution

        Parameter:
        -------
        None

        Return
        --------
        None
    '''
    # solve the ode numerically
    # Initialise paremeters for the ode model
    p0 = 1
    a1 = 2
    b = 1
    t0 = 0
    t1 = 20
    dt = 1/20
    nt = int(np.ceil((t1-t0)/dt))		# compute number of Euler steps to take
    t = t0+np.arange(nt+1)*dt
    q1 = 0*t
    q1.fill(1)
    q2 = np.sin(t)
    pars1 = [a1, b, p0]  # The ode model is dp/dt=-2(sin(t)-t)-(P-P0)
    ts, xsn = solve_ode_pressure(
        pressure_ode_model, t0, t1, dt, q1, q2, p0, pars1)
    # solve the ode with the given parameter analytically with wolfram alpha
    xsa = -3*np.exp(-t)-np.sin(t)+np.cos(t)+3
    # subplot the 3 graphs
    f, (ax1, ax2, ax3) = plt.subplots(1, 3)
    # First plot the graph for comparing analytical and numerical solution
    ax1.plot(ts, xsn, '-x')
    ax1.plot(ts, xsa, '-')
    ax1.title.set_text('Analytical solution matching')
    # Secondly plot the graph for absolute erro term
    ax2.plot(ts, np.absolute(xsn-xsa)/np.absolute(xsa), '-')
    ax2.title.set_text('Absolute Errors')
    # Finally plot convergence testing
    dtc = np.linspace(1, 7, 21)
    lastval = 0*dtc
    for i in range(21):
        # compute number of Euler steps to take
        ndt = 1/dtc[i]
        ntc = int(np.ceil((t1-t0)/ndt))
        tc = t0+np.arange(ntc+1)*ndt
        q1c = 0*tc
        q2c = np.sin(tc)
        tsc, xsnc = solve_ode_pressure(
            pressure_ode_model, t0, t1, ndt, q1c, q2c, p0, pars1)
        lastval[i] = xsnc[-1]
    ax3.plot(dtc, lastval, 'o')
    ax3.title.set_text('Convergence testing')
    plt.show()


def temperature_benchmark():
    ''' Check if the temperature fit the model

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
        tsc, xsnc = solve_ode(pressure_ode_model, t0=0, t1=10,
                              dt=1/dtc[i], x0=0, pars=[-1, 1, 1, 0])
        lastval[i] = xsnc[-1]
    ax3.plot(dtc, lastval, 'x')
    plt.show()


if __name__ == '__main__':
    pressure_benchmark()
