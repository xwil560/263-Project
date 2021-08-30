# ENGSCI 263: Lumped Parameter Model

import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.optimize import curve_fit

def pressure_ode_model(t, P, q1, q2, a, b, P0):
    ''' Return the derivative dP/dt at time, t(days), for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable, days
        P : float
            Dependent variable, kPa
        q1 : float
            Mass source rate of steam injection, tonnes/day
        q2 : float
            Mass sink rate of oil and water, m^3/day
        a : float
            Source/sink strength parameter.
        b : float
            Recharge strength parameter.
        P0 : float
            Ambient value of dependent variable, kPa.

        Returns:
        --------
        dPdt : float
            Derivative of dependent variable with respect to independent variable.

    '''
    #Pressure derivative calculation
    dPdt = -a*(3*q2-q1) - b*(P-P0)

    return dPdt

def temp_ode_model(t, T, q1, P, a, b, bt, P0, T0, M0):
    ''' Return the derivative dP/dt at time, t, for given parameters.

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
        P0 : float
            Ambient value of solved dependent variable, kPa.
        T0 : float
            Ambient value of dependent variable, DegreesC.
        M0 : float
            Initial mass in system, tonnes

        Returns:
        --------
        dTdt : float
            Derivative of dependent variable with respect to independent variable. degreesC/day

    '''

    #Checking direction of flow to determine temperature 
    if P > P0:
        Td = T 
    else:
        Td = T0

    #calculating the derivative for temperature
    dTdt = (q1/M0)*T - b/(a*M0)*(P-P0)*(Td-T) - bt*(T-T0)/M0

    return dTdt

def load_data():
    ''' read in and return the pilot data for the system

        Parameters:
        -----------
        None

        Returns:
        --------
        data : List
            list containing the pilot data for pressure, temperature, water and oil production rate, steam injection rate
            and the corresponding times for these data points.
    '''

    #Pressure data
    tP = np.genfromtxt('tr_p.txt',skip_header=1,usecols=0,delimiter=',')
    P = np.genfromtxt('tr_p.txt',skip_header=1,usecols=1,delimiter=',')
    #Temperature data
    tT = np.genfromtxt('tr_T.txt',skip_header=1,usecols=0,delimiter=',')
    T = np.genfromtxt('tr_T.txt',skip_header=1,usecols=1,delimiter=',')
    #water production data
    tW = np.genfromtxt('tr_water.txt',skip_header=1,usecols=0,delimiter=',')
    W = np.genfromtxt('tr_water.txt',skip_header=1,usecols=1,delimiter=',')
    #Oil production data
    tO = np.genfromtxt('tr_oil.txt',skip_header=1,usecols=0,delimiter=',')
    O = np.genfromtxt('tr_oil.txt',skip_header=1,usecols=1,delimiter=',')
    #Steam injection data
    tS = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=0,delimiter=',')
    S = np.genfromtxt('tr_steam.txt',skip_header=1,usecols=1,delimiter=',')

    #combining the data into a list to return
    data = [tP, P, tT, T, tW, W, tO, O, tS, S]

    return data

def interpolate_mass_source(t):
    ''' Return mass source parameter q1 for steam injection.

        Parameters:
        -----------
        t : array-like
            Vector of times (days) at which to interpolate the mass source.

        Returns:
        --------
        q1 : array-like
            Mass source of steam injection (tonnes per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    #navigating to data directory
    os.chdir("../data")

    #reading in observed data
    data = load_data()
    #selecting the injected steam rate and time 
    tS, S = data[8], data[9]
    #interpolating this data to find the steam injection at given times
    q1 = np.interp(t, tS, S)

    return q1

def interpolate_mass_sink(t):
    ''' Return mass sink parameter q2 for water and oil retrieval.

        Parameters:
        -----------
        t : array-like
            Vector of times(days) at which to interpolate the mass sink.

        Returns:
        --------
        q : array-like
            Mass sink rate of water and oil (m^3 per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    #changing working directory to data pathway
    os.chdir("../data")

    #using load_data function to read in all the experiment data
    data = load_data()
    #selecting the water and oil production data
    tW, W, tO, O = data[4], data[5], data[6], data[7]

    #interpolating water and oil production rates at the given times 
    water = np.interp(t, tW, W)
    oil = np.interp(t, tO, O)

    #combining the two sink rates into one total mass sink rate 
    q2 = water + oil

    return q2

def fit_pressure(t, a, b):
    '''function to solve for pressure (kPa) to use for curve fit parameter estimation

        Parameters:
        -----------
        t : array-like
            Vector of times(days) at which to interpolate the mass sink and source.
        a : float
            Source/sink strength parameter.
        b : float
            Recharge strength parameter.

        Returns:
        --------
        P : array-like
            pressure values of system at given times using these paramaters
        

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    
    #finding souce and sink values at given times 
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    #solving for pressure given a and b paramater values
    _,P = solve_ode_pressure(pressure_ode_model,0,216,1,q1,q2,1291.76,[a,b,1291.76])

    return P
    

def fit_temp(t, bt, M0):
    '''function to solve for Temperature (degC) to use for curve fit parameter estimation

        Parameters:
        -----------
        t : array-like
            Vector of times(days) at which to interpolate the mass sink and source.
        bt : float
            Temperature recharge strength parameter.
        M0 : float
            Initial mass in system, tonnes

        Returns:
        --------
        T : array-like
            Temperature of system at times t using these paramaters.

        Notes:
        ------
        Data interpolation can only be used between 0 - 216 days.

    '''
    #interpolating source and sink rates at times t
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    #a and b paramater values found using best fit for pressure
    a = 0.10414529 
    b = 0.05136862

    #solving for pressure 
    _,P = solve_ode_pressure(pressure_ode_model,0,216,1,q1,q2,1291.76,[a,b,1291.76])

    #Solving for temperature at given times, to return
    _,T = solve_ode_temp(temp_ode_model,0,216,1,q1,180.698,P,[a,b,bt,1291.76,180.698,M0])
    return T


def solve_ode_pressure(f, t0, t1, dt, q1, q2, P0, pars):
    ''' Solve the pressure ODE numerically.

        Parameters:
        -----------
        f : callable
            Function that returns dP/dt given variable and parameter inputs.
        t0 : float
            Initial time of solution, (day).
        t1 : float
            Final time of solution, (day).
        dt : float
            Time step length.
        P0 : float
            Initial value of solution, kPa.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        t : array-like
            Independent variable solution vector (days).
        P : array-like
            Dependent variable solution vector, kPa.

        Notes:
        ------
        ODE is solved using the Improved Euler Method.

    '''
    #the number of iterations needed
    iterations = int(np.ceil((t1-t0)/dt))

    #creating array for times of each iteration
    t = t0 + np.arange(iterations+1)*dt
    t[0] = t0
    for i in range(iterations):
        t[i+1] = t[i] + dt

    #creating array of system pressure at each time in t
    P = 0.*t
    P[0] = P0 #setting initial pressure value

    #solving for pressure using improved euler method
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
            Initial time of solution, days.
        t1 : float
            Final time of solution, days.
        dt : float
            Time step length, days.
        T0 : float
            Initial value of solution, degC.
        P : Array-like
            Array of values of previously solved P, kPa.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        t : array-like
            Independent variable solution vector, days.
        T : array-like
            Dependent variable solution vector, degC.

    '''
    #calculating number of interations required
    iterations = int(np.ceil((t1-t0)/dt))

    #creating array of times to solve for temperature at
    t = t0 + np.arange(iterations+1)*dt
    t[0] = t0
    for i in range(iterations):
        t[i+1] = t[i] + dt

    #creating empty array for Temerature at each time, with initial value
    T = 0.*t
    T[0] = T0

    #Solving for temperature on each day using improved euler method.
    for i in range(iterations):

        k1 = f(t[i], T[i], q1[i], P[i], *pars)
        k2 = f(t[i]+dt, T[i]+dt*k1, q1[i], P[i], *pars)
        T[i+1] = T[i] + 0.5*dt*(k1+k2)
        
    return t, T

def plot_models():
    ''' Plot the two pressure and temperature LPMs over top of the observed data.

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
    #changing working directory to data pathway
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

    #using curvefit to find optimum a and b values for pressure LPM
    Pf,_ = curve_fit(fit_pressure, t, Pi, [a,b])
    print(Pf)
    a = Pf[0]
    b = Pf[1]

    #using curvefit to find optimum b and initial mass values for temperautre LMP
    Tf,_ = curve_fit(fit_temp, t, Ti, [bt,M0])
    bt = Tf[0]
    M0 = Tf[1]
    
    # Initialising parameter arrays
    pars_P = [a, b, P0]
    pars_T = [a, b, bt, P0, T0, M0]

    #final solve for temperature and pressure over time using best fit paramaters
    tP, P = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1, q2, P0, pars_P)
    tT, T = solve_ode_temp(temp_ode_model, t0, t1, dt, q1, T0, P, pars_T)


    plt.rcParams["figure.figsize"] = (9, 6)
    f, ax1 = plt.subplots(1, 1) # Creating plot figure and axes
    ax2 = ax1.twinx() # Creating separate axis

    #plotting our LMPs and the given data
    ax1.plot(tP, P, 'k-', label='Pressure Best Fit')
    ax1.plot(tPe, Pe, 'k.', label='Data')
    ax2.plot(tT, T, 'r-', label='Temperature Best Fit')
    ax2.plot(tTe, Te, 'r.', label=' Data')

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
