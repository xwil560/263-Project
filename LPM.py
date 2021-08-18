# ENGSCI 263: Lumped Parameter Model

import numpy as np
from matplotlib import pyplot as plt
import os

def pressure_ode_model(t, P, q1, q2, a1, a2, b, P0):
    ''' Return the derivative dP/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable.
        P : float
            Dependent variable.
        q : float
            Source/sink rate.
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
    dPdt = (a1*q1 - a2*q2) - b*(P-P0)

    return dPdt

def temp_ode_model(t, T, q1, q2, m, P, a1, a2, b, d, P0, T0):
    ''' Return the derivative dP/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable.
        T : float
            Dependent variable.
        q : float
            Source/sink rate.
        m : float
            Mass variable
        P : float
            Pressure variable
        a : float
            Source/sink strength parameter.
        b : float
            Pressure recharge strength parameter.
        d : float
            Temperature recharge strength parameter.
        P0 : float
            Ambient value of solved dependent variable.
        T0 : float
            Ambient value of dependent variable.

        Returns:
        --------
        dTdt : float
            Derivative of dependent variable with respect to independent variable.

    '''
    # dTdt = (T*q/m - T*q2/m) - b/(a*m)*T*(P-P0) + d(T-T0)
    dTdt = (T*q1/m*a2 - T*q2/m*a1) - b/(a1*a2*m)*T*(P-P0) + d*(T-T0)
    
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
            Mass source (tonnes per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    os.chdir("../data")

    time_steam = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=0,delimiter=',')
    steam = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=1,delimiter=',')

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
            Mass sink (m^3 per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    os.chdir("../data")

    time_water = np.genfromtxt('tr_water.txt',skip_header=1,usecols=0,delimiter=',')
    water = np.genfromtxt('tr_water.txt',skip_header=1,usecols=1,delimiter=',')

    time_oil = np.genfromtxt('tr_oil.txt',skip_header=1,usecols=0,delimiter=',')
    oil = np.genfromtxt('tr_oil.txt',skip_header=1,usecols=1,delimiter=',')

    W = np.interp(t, time_water, water)
    O = np.interp(t, time_oil, oil)

    q2 = W + O

    return q2

def interpolate_mass_parameter(t):
    ''' Return mass parameter m for steam mass injected into reservoir.

        Parameters:
        -----------
        t : array-like
            Vector of times at which to interpolate the mass sink.

        Returns:
        --------
        m : array-like
            Steam mass (tonnes per day) interpolated at t.

        Notes:
        ------
        Data ends at 160 days, as steam injection stops, for the ease of interpolation
        I added a couple rows of zeros to take it to 216 days like the rest of the data.

        Data interpolation can only be used between 0 - 216 days.

    '''
    os.chdir("../data")

    time_steam = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=0,delimiter=',')
    steam = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=1,delimiter=',')

    m = np.interp(t, time_steam, steam)

    return m

def solve_ode_pressure(f, t0, t1, dt, P0, pars):
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
    for i in range (iterations):
        t[i+1] = t[i] + dt

    P = 0.*t
    P[0] = P0

    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

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
    for i in range (iterations):
        t[i+1] = t[i] + dt

    T = 0.*t
    T[0] = T0

    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)
    m = interpolate_mass_parameter(t)

    for i in range(iterations):

        k1 = f(t[i], T[i], q1[i], q2[i], m[i], P[i], *pars)
        k2 = f(t[i]+dt, T[i]+dt*k1, q1[i], q2[i], m[i], P[i], *pars)
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

    # Experimental data for Temperature and Pressure
    Pexp_t = np.genfromtxt('tr_p.txt',skip_header=1,usecols=0,delimiter=',')
    Pexp = np.genfromtxt('tr_p.txt',skip_header=1,usecols=1,delimiter=',')
    Texp_t = np.genfromtxt('tr_T.txt',skip_header=1,usecols=0,delimiter=',')
    Texp = np.genfromtxt('tr_T.txt',skip_header=1,usecols=1,delimiter=',')
    
    # PO = 381.187     1291.76
    P0 = 1291.76
    T0 = 180.698

    # Couldn't find much on specific values, well need to use SciPy curvefit.
    # 9.8 * 10**-3
    # a1 should be smaller than a2, retrieval of fluids need to be weighted heavier.
    a1 = 2
    a2 = 3.75
    b = 1
    d = 100

    pars_P = [a1, a2, b, P0]
    pars_T = [a1, a2, b, d, P0, T0]

    P_t, P = solve_ode_pressure(pressure_ode_model, 0, 216, 1, 1291.76, pars_P)
    #T_t, T = solve_ode_temp(temp_ode_model, 0, 216, 1, T0, P, pars_T)

    f, ax1 = plt.subplots(1,1) # Creating plot figure and axes
    ax2 = ax1.twinx() # Creating separate axis

    ax1.plot(Pexp_t, Pexp, 'r.', label='EXP PRESSURE')
    ax1.plot(P_t, P, 'r-', label='INTERPOLATED P')
    #ax2.plot(Texp_t, Texp, 'b.', label='EXP TEMP')
    #ax2.plot(T_t, T, 'b-', label='INTERPOLATED T')

    # Setting y limits for each axes, drawing labels and legends 
    ax1.set_ylabel('Pressure')
    ax2.set_ylabel('Temperature')
    ax1.set_xlabel('Time (days)')
    ax1.set_title('Pressure and Temperature Models')

    # Setting colour of axes
    ax1.tick_params(axis='y',colors='red')
    ax1.yaxis.label.set_color('red')
    ax2.tick_params(axis='y',colors='blue')
    ax2.yaxis.label.set_color('blue')

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
    


