# ENGSCI 263: Lumped Parameter Model

import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.optimize import curve_fit

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
    dPdt = -a*(3*q2-q1) - b*(P-P0)

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
        Td = T
    else:
        Td = T0

    dTdt = (q1/M0)*T - b/(a*M0)*(P-P0)*(Td-T) - bt*(T-T0)/M0

    return dTdt

def load_data():

    tP = np.genfromtxt('tr_p.txt',skip_header=1,usecols=0,delimiter=',')
    P = np.genfromtxt('tr_p.txt',skip_header=1,usecols=1,delimiter=',')

    tT = np.genfromtxt('tr_T.txt',skip_header=1,usecols=0,delimiter=',')
    T = np.genfromtxt('tr_T.txt',skip_header=1,usecols=1,delimiter=',')

    tW = np.genfromtxt('tr_water.txt',skip_header=1,usecols=0,delimiter=',')
    W = np.genfromtxt('tr_water.txt',skip_header=1,usecols=1,delimiter=',')

    tO = np.genfromtxt('tr_oil.txt',skip_header=1,usecols=0,delimiter=',')
    O = np.genfromtxt('tr_oil.txt',skip_header=1,usecols=1,delimiter=',')

    tS = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=0,delimiter=',')
    S = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=1,delimiter=',')

    data = [tP, P, tT, T, tW, W, tO, O, tS, S]

    return data

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

    data = load_data()
    tS, S = data[8], data[9]

    q1 = np.interp(t, tS, S)

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

    data = load_data()
    tW, W, tO, O = data[4], data[5], data[6], data[7]

    water = np.interp(t, tW, W)
    oil = np.interp(t, tO, O)

    q2 = water + oil

    return q2

def fit_pressure(t, a, b):
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    _,P = solve_ode_pressure(pressure_ode_model,0,216,1,q1,q2,1291.76,[a,b,1291.76])
    return P
    
def fit_temp(t, bt, M0):
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    a = 0.10414529 
    b = 0.05136862
    
    _,P = solve_ode_pressure(pressure_ode_model,0,216,1,q1,q2,1291.76,[a,b,1291.76])

    _,T = solve_ode_temp(temp_ode_model,0,216,1,q1,180.698,P,[a,b,bt,1291.76,180.698,M0])
    return T

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

def solve_ode_temp(f, t0, t1, dt, q1, T0, P, pars):
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

    # Experimental data for pressure and temperature
    data = load_data()
    tPe, Pe, tTe, Te = data[0], data[1], data[2], data[3]

    # Initialising time array
    t0 = 0
    t1 = 216
    dt = 1
    
    iterations = int(np.ceil((t1-t0)/dt))
    t = t0 + np.arange(iterations+1)*dt

    # Initial values of pressure and temperature
    P0 = Pe[0]
    T0 = Te[0]
    
    # Initial guesses for parameters
    a = 0.2
    b = 0.1
    bt = 0.8
    M0 = 50000

    # Calling the interpolation functions for q1 and q2 arrays
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    # Fit parameters
    Pi = np.interp(t, tPe, Pe)
    Ti = np.interp(t, tTe, Te)

    Pf,_ = curve_fit(fit_pressure, t, Pi, [a,b])
    print(Pf)
    a = Pf[0]
    b = Pf[1]

    Tf,_ = curve_fit(fit_temp, t, Ti, [bt,M0])
    bt = Tf[0]
    M0 = Tf[1]
    
    # Initialising parameter arrays
    pars_P = [a, b, P0]
    pars_T = [a, b, bt, P0, T0, M0]

    tP, P = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1, q2, P0, pars_P)
    tT, T = solve_ode_temp(temp_ode_model, t0, t1, dt, q1, T0, P, pars_T)

    plt.rcParams["figure.figsize"] = (9, 6)
    f, ax1 = plt.subplots(1, 1) # Creating plot figure and axes
    ax2 = ax1.twinx() # Creating separate axis

    ax1.plot(tP, P, 'k-', label='Pressure Best Fit')
    ax1.plot(tPe, Pe, 'k.', label='Data')
    ax2.plot(tT, T, 'r-', label='Temperature Best Fit')
    ax2.plot(tTe, Te, 'r.', label=' Data')

    # Setting y limits for each axes, drawing labels and legends
    ax1.set_ylabel('Pressure (Pa)')
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
