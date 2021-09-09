# imports
import numpy as np
from matplotlib import pyplot as plt
from main import *


def pressure_benchmark():
    ''' Check if the pressure ode is implemented correctly by comparing it to the analytical solution

        Parameter:
        -------
        None

        Return
        --------
        None
    '''

    # Initialise paremeters for the ode model
    p0 = 1
    a1 = 2
    b = 1
    t0 = 0
    t1 = 20
    dt = 1/20

    nt = int(np.ceil((t1-t0)/dt))  # Compute number of Euler steps to take
    t = t0+np.arange(nt+1)*dt
    q1 = 0*t
    q1.fill(1)
    q2 = np.sin(t)
    pars1 = [a1, b, p0]  # The ode model is dp/dt=-2(sin(t)-t)-(P-P0)
    # solve the ode numerically
    ts, xsn = solve_ode_pressure(
        pressure_ode_model, t0, t1, dt, q1, q2, p0, pars1)
    # Solve the ode with the given parameter analytically with Wolfram Alpha
    xsa = -3*np.exp(-t)-np.sin(t)+np.cos(t)+3

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


def temp_model(t, T, q1, q2, P, a, b, bt, Pa, Ta, M0):
    ''' A new temperature without the if commands to solve it analytically
    Parameters:
        -----------
        t : float
            Independent variable, days.
        T : float
            Dependent variable, degreesC
        q1 : float
            Mass source rate of steam injection, tonnes/day
        P : float
            Pressure variable, kPa.
        a : float
            Source/sink strength parameter.
        b : float
            Pressure recharge strength parameter.
        bt : float
            Temperature recharge strength parameter.
        Pa : float
            Ambient value of solved dependent variable, kPa.
        Ta : float
            Ambient value of dependent variable, DegreesC.
        M0 : float
            Initial mass in system, tonnes

        Returns:
        --------
        dTdt : float
            Derivative of dependent variable with respect to independent variable. degreesC/day

    '''

    dTdt = ((q1-q2)/M0)*(500-T) - (b/(a*M0))*(P-Pa)*T - bt*(T-Ta)

    return dTdt


def temperature_benchmark():
    ''' Check if the temperature fit the model

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
    dt = 1
    bt = 0.5
    T0 = 25
    M0 = 100

    nt = int(np.ceil((t1-t0)/dt))  # Compute number of Euler steps to take
    t = t0+np.arange(nt+1)*dt
    q1 = 0*t
    q1.fill(1)
    q2 = 11*q1
    P = 150*q1
    # The ode model is 0.15T-37.5
    pars2 = [a1, b, bt, p0, T0, M0]
    # solve the ode numerically
    ts, xsn = solve_ode_temp(temp_model, t0, t1, dt, q1, T0, P, pars2)
    # Solve the ode with the given parameter analytically with Wolfram Alpha
    xsa = 275*np.exp(-0.15*t)-250

    # Subplot the 3 graphs
    f, (ax1, ax2, ax3) = plt.subplots(1, 3)

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
        q1c.fill(1)
        q2c = 11*q1c
        Pc = 150*q1c
        tsc, xsnc = solve_ode_temp(
            temp_model, t0, t1, ndt, q1c, q2c, T0, Pc, pars2)
        lastval[i] = xsnc[-1]

    ax3.plot(dtc, lastval, 'o')
    ax3.title.set_text('Timestep Convergence')
    ax3.set_ylabel("X(t)")
    ax3.set_xlabel("1/Δt")

    plt.tight_layout()
    plt.show()


if __name__ == '__main__':
    pressure_benchmark()
