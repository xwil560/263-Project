# ENGSCI 263: Lumped Parameter Model

import numpy as np
from matplotlib import pyplot as plt

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
    dPdt = (a1*q1) - (a2*q2) - b*(P-P0)

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
    dTdt = (T*q1/(m*a2) - T*q2/(m*a1)) - b/(a1*a2*m)*T*(P-P0) + d(T-T0)
    
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
    time_steam = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=0,delimiter=',')
    steam = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=1,delimiter=',')

    m = np.interp(t, time_steam, steam)

    return m

def solve_ode_pressure(f, t0, t1, dt, P0, pars):
    ''' Solve the pressure ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dPdt given variable and parameter inputs.
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
    # time array
    t = range(t0,t1,dt)
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)
    # amount of euler steps to complete
    n = len(t)
    # initialise solution array
    P = np.zeros(n)
    P[0] = P0

    i = 0
    while (i < n-1):

        fn1 = f(t[i], P[i], q1[i], q2[i], *pars)
        Pk1 = P[i] + dt*fn1
        fn2 = f(t[i+1], Pk1, q1[i], q2[i], *pars)
        P[i+1] = P[i] + dt*(fn1 + fn2)/2
        
        i = i+1

    return t, P

def solve_ode_temp(f, t0, t1, dt, T0, P, pars):
    ''' Solve the temperature ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dxdt given variable and parameter inputs.
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
    # time array
    t = range(t0,t1,dt)
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)
    m = interpolate_mass_parameter(t) # need to fix this???
    '''
    might have to make all values of m that equal 0 to be 1 instead, because of division by zero. 
    '''
    i = 0
    while (i < len(t)):
        if m[i] < 1:
            m[i] = 1
        i = i+1

    # amount of euler steps to complete
    n = len(t)
    # initialise solution array
    T = np.zeros(n)
    T[0] = T0
    i = 0
    
    while (i < n-1):

        fn1 = f(t[i], T[i], q1[i], q2[i], m[i], P[i], *pars)
        xk1 = T[i] + dt*fn1
        fn2 = f(t[i+1], xk1, q1[i], q2[i], m[i], P[i] *pars)
        T[i+1] = T[i] + dt*(fn1 + fn2)/2
        
        i = i+1

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
    # Experimental data for Temperature and Pressure
    texp = np.genfromtxt('tr_T.txt',skip_header=1,usecols=0,delimiter=',')
    Texp = np.genfromtxt('tr_T.txt',skip_header=1,usecols=1,delimiter=',')

    texp2 = np.genfromtxt('tr_p.txt',skip_header=1,usecols=0,delimiter=',')
    Pexp = np.genfromtxt('tr_p.txt',skip_header=1,usecols=1,delimiter=',')
    
    # pO = 381.187     1291.76
    p0 = 1291.76
    T0 = 180.698

    # Couldn't find much on specific values, well need to use SciPy curvefit.
    # 9.8 * 10**-3
    # a1 should be smaller than a2, retrieval of fluids need to be weighted heavier.
    a1 = 2
    a2 = 3.75
    b = 1
    d = 100

    pars1 = [a1, a2, b, p0]
    pars2 = [a1, a2, b, d, p0, T0]

    t1, P = solve_ode_pressure(pressure_ode_model, 0, 216, 1, 1291.76, pars1)
    #t2, T = solve_ode_temp(temp_ode_model, 1, 216, 1, T0, P, pars2)

    f,ax1 = plt.subplots(nrows=1,ncols=1) # creates plot figure and axes

    ax1.plot(texp2, Pexp, 'r.', label='EXP PRESSURE')
    #ax1.plot(texp, Texp, 'g.', label='EXP TEMP')
    ax1.plot(t1, P, 'b-', label='INTERPOLATED P')
    #ax1.plot(t2, T, 'b-', label='INTERPOLATED T')

    # Setting y limits for each axes, drawing labels and legends 
    ax1.set_ylabel('Pressure')
    ax1.set_title('Pressure and Temperature Models')

	# Either show the plot to the screen or save a version of it to the disk
    save_figure = False
    
    if not save_figure:
        plt.show()
    else:
        plt.savefig('Pressure_Temp.png',dpi=300)

    return

if __name__ == "__main__":
    plot_models()
    


