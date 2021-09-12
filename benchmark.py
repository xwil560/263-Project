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

    # Initialise parameters for the ode model
    p0 = 1.75
    a1 = -2
    b = 1/1.75
    t0 = 0
    t1 = 20
    dt = 1/20

    nt = int(np.ceil((t1-t0)/dt))  # Compute number of Euler steps to take
    t = t0+np.arange(nt+1)*dt
    q1 = 0*t
    q1.fill(1)
    q2 = np.sin(t)
    pars1 = [a1, b, p0]  # The ode model is dp/dt = P + 2sin(t) - 1
    # Solve the ode numerically
    ts, xsn = solve_ode_pressure(
        pressure_ode_model, t0, t1, dt, q1, q2, p0, pars1)
    # Solve the ode with the given parameter analytically with Wolfram Alpha
    xsa = 3.75*np.exp(-t)+np.sin(t)-np.cos(t)-1

    # Creating plot figure and axes
    plt.rcParams["figure.figsize"] = (15, 5)
    _, (ax1, ax2, ax3) = plt.subplots(1, 3)

    # Plotting the graph for comparing analytical and numerical solution
    ax1.plot(ts, xsn, '-x')
    ax1.plot(ts, xsa, '-')
    ax1.title.set_text('Pressure Benchmark')
    ax1.set_ylabel("Pressure (kPa)")
    ax1.set_xlabel("Time (days)")

    # Plotting the graph for absolute error term
    ax2.plot(ts, np.absolute(xsn-xsa)/np.absolute(xsa), '-')
    ax2.title.set_text('Error Analysis')
    ax2.set_ylabel("Relative Error Against Benchmark")
    ax2.set_xlabel("Time (days)")

    # Plotting the graph for convergence testing
    dtc = np.linspace(1, 7, 21)
    lastval = 0*dtc

    for i in range(21):
        # Compute number of Euler steps to take
        ndt = 1/dtc[i]
        ntc = int(np.ceil((t1-t0)/ndt))
        tc = t0+np.arange(ntc+1)*ndt
        q1c = 0*tc
        q2c = np.sin(tc)
        tsc, xsnc = solve_ode_pressure(
            pressure_ode_model, t0, t1, ndt, q1c, q2c, p0, pars1)
        lastval[i] = xsnc[-1]

    ax3.plot(dtc, lastval, 'o')
    ax3.title.set_text('Timestep Convergence')
    ax3.set_ylabel("X(t)")
    ax3.set_xlabel("1/Δt")

    plt.tight_layout()
    plt.show()


def temperature_benchmark():
    ''' Check if the temperature fit the model in the first scenario of the solver

        Parameter:
        ---------
        None

        Return:
        -------
        None
    '''
    # Initialise paremeters for the ode model
    p0 = 200
    a1 = 2
    b = 1
    t0 = 0
    t1 = 20
    dt = 1/4
    bt = 0.5
    T0 = 25
    M0 = 100

    nt = int(np.ceil((t1-t0)/dt))  # Compute number of Euler steps to take
    t = t0+np.arange(nt+1)*dt
    q1 = 0*t
    q1.fill(10)
    P = 25*q1
    # The ode model is dT/dt = 65.815 - 0.1T
    pars2 = [a1, b, bt, p0, T0, M0]
    # solve the ode numerically
    ts, xsn = solve_ode_temp(temp_ode_model, t0, t1, dt, q1, T0, P, pars2)
    # Solve the ode with the given parameter analytically with Wolfram Alpha
    xsa = 658.15-633.15*np.exp(-0.1*t)

    # Creating plot figure and axes
    plt.rcParams["figure.figsize"] = (15, 5)
    _, (ax1, ax2, ax3) = plt.subplots(1, 3)

    # Plotting the graph for comparing analytical and numerical solution
    ax1.plot(ts, xsn, '-x')
    ax1.plot(ts, xsa, '-')
    ax1.title.set_text('Temperature Benchmark')
    ax1.set_ylabel("Temperature(degrees C)")
    ax1.set_xlabel("Time (days)")

    # Plotting the graph for absolute error term
    ax2.plot(ts, np.absolute(xsn-xsa)/np.absolute(xsa), '-')
    ax2.title.set_text('Error Analysis')
    ax2.set_ylabel("Relative Error Against Benchmark")
    ax2.set_xlabel("Time (days)")

    # Plotting the graph for convergence testing
    dtc = np.linspace(1, 20, 21)
    lastval = 0*dtc
    for i in range(21):
        # Compute number of Euler steps to take
        ndt = 1/dtc[i]
        ntc = int(np.ceil((t1-t0)/ndt))
        tc = t0+np.arange(ntc+1)*ndt
        q1c = 0*tc
        q1c.fill(10)
        Pc = 25*q1c
        tsc, xsnc = solve_ode_temp(
            temp_ode_model, t0, t1, ndt, q1c, T0, Pc, pars2)
        lastval[i] = xsnc[-1]

    ax3.plot(dtc, lastval, 'o')
    ax3.title.set_text('Timestep Convergence')
    ax3.set_ylabel("X(t)")
    ax3.set_xlabel("1/Δt")

    plt.tight_layout()
    plt.show()
