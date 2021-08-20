# ENGSCI 263: Lumped Parameter Model

import numpy as np
from matplotlib import pyplot as plt
import os


def pressure_ode_model(t, P, q1, q2, a, b, P0):
    ''' Return the derivative dP/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable.
        P : float
            Dependent variable.
        q1 : float
            Mass source rate of steam injection.
        q2 : float
            Mass sink rate of oil and water.
        a : float
            Source/sink strength parameter.
        b : float
            Recharge strength parameter.
        P0 : float
            Ambient value of dependent variable.

        Returns:
        --------
        dPdt : float
            Derivative of dependent variable with respect to independent variable.

    '''
    dPdt = -a*(q2-q1) - b*(P-P0)

    return dPdt


def temp_ode_model(t, T, q1, P, a, b, bt, P0, T0, M0):
    ''' Return the derivative dP/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable.
        T : float
            Dependent variable.
        q1 : float
            Mass source rate of steam injection.
        P : float
            Pressure variable.
        a : float
            Source/sink strength parameter.
        b : float
            Pressure recharge strength parameter.
        bt : float
            Temperature recharge strength parameter.
        P0 : float
            Ambient value of solved dependent variable.
        T0 : float
            Ambient value of dependent variable.
        M0 : float
            Initial mass in system.

        Returns:
        --------
        dTdt : float
            Derivative of dependent variable with respect to independent variable.

    '''
    if P > P0:
        td = T
    else:
        td = T0

    dTdt = (q1/M0)*T - b/(a*M0) * (P-P0) * (td-T) + bt*(T-T0)

    return dTdt


def interpolate_mass_source(t):
    ''' Return mass source parameter q1 for steam injection.

        Parameters:
        -----------
        t : array-like
            Vector of times at which to interpolate the mass source.

        Returns:
        --------
        q1 : array-like
            Mass source of steam injection (tonnes per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    os.chdir("../data")

    time_steam = np.genfromtxt(
        'tr_steam.txt', skip_header=1, usecols=0, delimiter=',')
    steam = np.genfromtxt('tr_steam.txt', skip_header=1,
                          usecols=1, delimiter=',')

    q1 = np.interp(t, time_steam, steam)

    return q1


def interpolate_mass_sink(t):
    ''' Return mass sink parameter q2 for water and oil retrieval.

        Parameters:
        -----------
        t : array-like
            Vector of times at which to interpolate the mass sink.

        Returns:
        --------
        q : array-like
            Mass sink rate of water and oil (m^3 per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    os.chdir("../data")

    time_water = np.genfromtxt(
        'tr_water.txt', skip_header=1, usecols=0, delimiter=',')
    water = np.genfromtxt('tr_water.txt', skip_header=1,
                          usecols=1, delimiter=',')

    time_oil = np.genfromtxt(
        'tr_oil.txt', skip_header=1, usecols=0, delimiter=',')
    oil = np.genfromtxt('tr_oil.txt', skip_header=1, usecols=1, delimiter=',')

    W = np.interp(t, time_water, water)
    O = np.interp(t, time_oil, oil)

    q2 = W + O

    return q2


def solve_ode_pressure(f, t0, t1, dt, q1, q2, P0, pars):
    ''' Solve the pressure ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dP/dt given variable and parameter inputs.
        t0 : float
            Initial time of solution.
        t1 : float
            Final time of solution.
        dt : float
            Time step length.
        P0 : float
            Initial value of solution.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        t : array-like
            Independent variable solution vector.
        P : array-like
            Dependent variable solution vector.

        Notes:
        ------
        ODE is solved using the Improved Euler Method.

    '''
    iterations = int(np.ceil((t1-t0)/dt))

    t = t0 + np.arange(iterations+1)*dt
    t[0] = t0
    for i in range(iterations):
        t[i+1] = t[i] + dt

    P = 0.*t
    P[0] = P0

    for i in range(iterations):

        k1 = f(t[i], P[i], q1[i], q2[i], *pars)
        k2 = f(t[i]+dt, P[i]+dt*k1, q1[i], q2[i], *pars)
        P[i+1] = P[i] + 0.5*dt*(k1+k2)

    return t, P


def solve_ode_temp(f, t0, t1, dt, T0, P, pars):
    ''' Solve the temperature ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dT/dt given variable and parameter inputs.
        t0 : float
            Initial time of solution.
        t1 : float
            Final time of solution.
        dt : float
            Time step length.
        T0 : float
            Initial value of solution.
        P : Array-like
            Array of values of previously solved P.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        t : array-like
            Independent variable solution vector.
        T : array-like
            Dependent variable solution vector.

    '''
    iterations = int(np.ceil((t1-t0)/dt))

    t = t0 + np.arange(iterations+1)*dt
    t[0] = t0
    for i in range(iterations):
        t[i+1] = t[i] + dt

    T = 0.*t
    T[0] = T0

    q1 = interpolate_mass_source(t)

    for i in range(iterations):

        k1 = f(t[i], T[i], q1[i], P[i], *pars)
        k2 = f(t[i]+dt, T[i]+dt*k1, q1[i], P[i], *pars)
        T[i+1] = T[i] + 0.5*dt*(k1+k2)
        
    return t, T


def plot_models():
    ''' Plot the kettle LPM over top of the data.

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":

        It contains commands to read and plot the experimental data, run and 
        plot the LPM for hard coded parameters, and then either display the 
        plot to the screen or save it to the disk.

    '''
    os.chdir("data")
    # Get the interpolation function to for the ode solver
    t0 = 0
    t1 = 216
    dt = 1
    t = range(t0, t1, dt)
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    # Experimental data for pressure and temperature
    Pt_exp = np.genfromtxt('tr_p.txt',skip_header=1,usecols=0,delimiter=',')
    P_exp = np.genfromtxt('tr_p.txt',skip_header=1,usecols=1,delimiter=',')
    Tt_exp = np.genfromtxt('tr_T.txt',skip_header=1,usecols=0,delimiter=',')
    T_exp = np.genfromtxt('tr_T.txt',skip_header=1,usecols=1,delimiter=',')
    
    # Initial values of pressure and temperature
    P0 = 1291.76
    T0 = 180.698
    
    # Initialising parameters
    a = 2.5
    b = 1
    bt = 0.00001
    M0 = 1000000
    
    # Initialising parameter arrays
    pars_P = [a, b, P0]
    pars_T = [a, b, bt, P0, T0, M0]

    Pt, P = solve_ode_pressure(
        pressure_ode_model, t0, t1, dt, q1, q2, P0, pars_P)
    Tt, T = solve_ode_temp(temp_ode_model, 0, 216, 1, T0, P, pars_T)

    plt.rcParams["figure.figsize"] = (9, 6)
    f, ax1 = plt.subplots(1, 1)  # Creating plot figure and axes
    ax2 = ax1.twinx()  # Creating separate axis

    ax1.plot(Pt, P, 'k-', label='Pressure Best Fit')
    ax1.plot(Pt_exp, P_exp, 'k.', label='Data')
    ax2.plot(Tt, T, 'r-', label='Temperature Best Fit')
    ax2.plot(Tt_exp, T_exp, 'r.', label=' Data')

    # Setting y limits for each axes, drawing labels and legends
    ax1.set_ylabel('Pressure (kPa)')
    ax2.set_ylabel('Temperature ($^{0}C$)')
    ax1.set_xlabel('Time (days)')
    ax1.set_title('Pressure and Temperature Models')
    ax1.legend(loc='upper right')
    ax2.legend(loc='lower left')

    # Either show the plot to the screen or save a version of it to the disk
    os.chdir("../plots")
    save_figure = False

    if not save_figure:
        plt.show()
    else:
        plt.savefig('Pressure_Temp.png', dpi=300)

    return


if __name__ == "__main__":
    plot_models()
