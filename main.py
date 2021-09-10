# ENGSCI 263: Lumped Parameter Model

# Imports
import numpy as np
from matplotlib import pyplot as plt
import os
from scipy.optimize import curve_fit


def pressure_ode_model(t, P, q1, q2, a, b, Pa):
    ''' Return the derivative dP/dt at time, t(days), for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable, days
        P : float
            Dependent variable, Pa
        q1 : float
            Mass source rate of steam injection, kg/day
        q2 : float
            Mass sink rate of oil and water, kg/day
        a : float
            Source/sink strength parameter.
        b : float
            Recharge strength parameter.
        Pa : float
            Ambient value of dependent variable, Pa.

        Returns:
        --------
        dPdt : float
            Derivative of dependent variable with respect to independent variable.

    '''
    # Calculating the derivative for pressure
    dPdt = -a*(q2-q1) - b*(1.75*P-Pa)

    return dPdt


def temp_ode_model(t, T, q1, P, a, b, bt, Pa, Ta, M0):
    ''' Return the derivative dP/dt at time, t, for given parameters.

        Parameters:
        -----------
        t : float
            Independent variable, days.
        T : float
            Dependent variable, degreesC
        q1 : float
            Mass source rate of steam injection, kg/day
        P : float
            Pressure variable, Pa.
        a : float
            Source/sink strength parameter.
        b : float
            Pressure recharge strength parameter.
        bt : float
            Temperature recharge strength parameter.
        Pa : float
            Ambient value of solved dependent variable, Pa.
        Ta : float
            Ambient value of dependent variable, Kelvins.
        M0 : float
            Initial mass in system, kgs.

        Returns:
        --------
        dTdt : float
            Derivative of dependent variable with respect to independent variable, Kelvins/day

    '''
    # Checking direction of flow to determine temperature
    if P > Pa:
        Tp = T
    else:
        Tp = Ta

    # Calculating the derivative for temperature
    dTdt = (q1/M0)*(533.15-T) - (b/(a*M0))*(P-Pa)*(Tp-3*T) - bt*(T-Ta)

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
    # Pressure data
    tP = np.genfromtxt('tr_p.txt', skip_header=1, usecols=0, delimiter=',')
    P = np.genfromtxt('tr_p.txt', skip_header=1, usecols=1, delimiter=',')
    # Temperature data
    tT = np.genfromtxt('tr_T.txt', skip_header=1, usecols=0, delimiter=',')
    T = np.genfromtxt('tr_T.txt', skip_header=1, usecols=1, delimiter=',')
    # Water production data
    tW = np.genfromtxt('tr_water.txt', skip_header=1, usecols=0, delimiter=',')
    W = np.genfromtxt('tr_water.txt', skip_header=1, usecols=1, delimiter=',')
    # Oil production data
    tO = np.genfromtxt('tr_oil.txt', skip_header=1, usecols=0, delimiter=',')
    O = np.genfromtxt('tr_oil.txt', skip_header=1, usecols=1, delimiter=',')
    # Steam injection data
    tS = np.genfromtxt('tr_steam.txt', skip_header=1, usecols=0, delimiter=',')
    S = np.genfromtxt('tr_steam.txt', skip_header=1, usecols=1, delimiter=',')

    # Combining the data into a list to return and converting into SI units
    data = [tP, P*1000, tT, T+273.15, tW, W*1000, tO, O*1000, tS, S*1000]

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
            Mass source of steam injection (kg per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 221 days.

    '''
    # Navigating to data directory
    os.chdir("../data")

    # Reading in observed data
    data = load_data()
    # Selecting the injected steam rate and time
    tS, S = data[8], data[9]

    # Interpolating this data to find the steam injection at given times
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
            Mass sink rate of water and oil (kg per day) interpolated at t.

        Notes:
        ------
        Data interpolation can only be used between 0 - 221 days.

    '''
    # Changing working directory to data pathway
    os.chdir("../data")

    # Using load_data function to read in all the experiment data
    data = load_data()
    # Selecting the water and oil production data
    tW, W, tO, O = data[4], data[5], data[6], data[7]

    # Interpolating water and oil production rates at the given times 
    water = np.interp(t, tW, W)
    oil = np.interp(t, tO, O)

    # Combining the two sink rates into one total mass sink rate
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
        Data interpolation can only be used between 0 - 221 days.

    '''
    # Finding souce and sink values at given times
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    # Solving for pressure given a and b paramater values
    _, P = solve_ode_pressure(pressure_ode_model, 0, 221, 1, q1, q2, 1291760, [a,b,1291760])

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
            Initial mass in system, kg

        Returns:
        --------
        T : array-like
            Temperature of system at times t using these paramaters.

        Notes:
        ------
        Data interpolation can only be used between 0 - 221 days.

    '''
    # Interpolating source and sink rates at times t
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    # Paramater values found using best fit for pressure
    a = 0.16229278 
    b = 0.02694718

    # Solving for pressure
    _, P = solve_ode_pressure(pressure_ode_model, 0, 221, 1, q1, q2, 1291760, [a,b,1291760])

    # Solving for temperature at given times to return
    _, T = solve_ode_temp(temp_ode_model, 0, 221, 1, q1, 453.848, P, [a,b,bt,1291760,453.848,M0])

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
        q1: array-like
            Steam injection mass rates (kg per day)
        q2: array-like
            Water and oil extraction rates (kg/day)
        P0 : float
            Initial value of solution, Pa.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        t : array-like
            Independent variable solution vector (days).
        P : array-like
            Dependent variable solution vector, Pa.

        Notes:
        ------
        ODE is solved using the Improved Euler Method.

    '''
    # Number of iterations needed
    iterations = int(np.ceil((t1-t0)/dt))

    # Creating array for times of each iteration
    t = t0 + np.arange(iterations+1)*dt
    t[0] = t0
    for i in range(iterations):
        t[i+1] = t[i] + dt

    # Creating array of system pressure at each time in t
    P = 0.*t
    P[0] = P0  # Setting initial pressure value

    # Solving for pressure using improved euler method
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
        q1: array-like
            Steam injection mass rates (kg per day)
        T0 : float
            Initial value of solution, Kelvins.
        P : Array-like
            Array of values of previously solved P, Pa.
        pars : array-like
            List of parameters passed to ODE function f.

        Returns:
        --------
        t : array-like
            Independent variable solution vector, days.
        T : array-like
            Dependent variable solution vector, Kelvins.

    '''    
    # Calculating number of interations required
    iterations = int(np.ceil((t1-t0)/dt))

    # Creating array of times to solve for temperature at
    t = t0 + np.arange(iterations+1)*dt
    t[0] = t0
    for i in range(iterations):
        t[i+1] = t[i] + dt

    # Creating empty array for Temperature at each time, with initial value
    T = 0.*t
    T[0] = T0

    # Solving for temperature on each day using improved euler method.
    for i in range(iterations):

        k1 = f(t[i], T[i], q1[i], P[i], *pars)
        k2 = f(t[i]+dt, T[i]+dt*k1, q1[i], P[i], *pars)
        T[i+1] = T[i] + 0.5*dt*(k1+k2)

    return t, T


def formulate_models():
    ''' Formulate the two pressure and temperature LPMs from the observed data.

        Parameters:
        -----------
        none

        Returns:
        --------
        t : ndarray
            times at which data is measured/interpolated, days
        tP : ndarray
            times at which pressure was solved at, days
        P : ndarray
            Pressure values solved at tP times, Pa
        tPe : ndarray
            times at which experiment pressure data was collected at, days
        Pe : ndarray
            Experiment Pressure values, Pa
        tT : ndarray
            times as which temperature was solved at, days
        T : ndarray
            Temperature values solved at tT times, Kelvins
        tTe : ndarray
            times at which experiment pressuer data was collected at, days
        Te : ndarray
            experiment Temerature values, Kelvins

    '''
    # Changing working directory to data pathway
    os.chdir("data")

    # Loading experimental data for pressure and temperature
    data = load_data()
    tPe, Pe, tTe, Te = data[0], data[1], data[2], data[3]

    # Initialising time array
    t0 = 0
    t1 = 221
    dt = 1
    iterations = int(np.ceil((t1-t0)/dt))
    t = t0 + np.arange(iterations+1)*dt

    # Initial values of pressure and temperature
    P0 = Pe[0]
    T0 = Te[0]

    #ambient pressure and temperature values 
    Pa = 1291760
    Ta = 453.848

    # Initial guesses for parameters
    a = 0.15
    b = 0.02
    bt = 0.9
    M0 = 4000000

    # Calling the interpolation functions for q1 and q2 arrays
    q1 = interpolate_mass_source(t)
    q2 = interpolate_mass_sink(t)

    # Fit parameters
    Pi = np.interp(t, tPe, Pe)
    Ti = np.interp(t, tTe, Te)

    # Using curvefit to find optimum a and b values for pressure LPM
    Pf,_ = curve_fit(fit_pressure, t, Pi, [a,b])
    a = Pf[0]
    b = Pf[1]
    
    # Using curvefit to find optimum b and initial mass values for temperautre LMP
    Tf,_ = curve_fit(fit_temp, t, Ti, [bt,M0])
    bt = Tf[0]
    M0 = Tf[1]

    print(Pf, Tf)

    # Initialising parameter arrays
    pars_P = [a, b, Pa]
    pars_T = [a, b, bt, Pa, Ta, M0]

    # Final solve for temperature and pressure over time using best fit paramaters
    tP, P = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1, q2, P0, pars_P)
    tT, T = solve_ode_temp(temp_ode_model, t0, t1, dt, q1, T0, P, pars_T)

    # Calculating misfits
    Pmf = np.full(len(t),0)
    Tmf = np.full(len(t),0)

    for i in range(len(t)):
        Pmf[i] = P[i] - Pi[i]
        Tmf[i] = T[i] - Ti[i]
        
    return t, tP, P, tPe, Pe, tT, T, tTe, Te, Pmf, Tmf


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

    # Loading data from model formulation function
    data = formulate_models()
    tP, P, tPe, Pe, tT, T, tTe, Te, Pmf, Tmf = data[1], data[2], data[3], data[4], data[5], data[6], data[7], data[8], data[9], data[10]

    plt.rcParams["figure.figsize"] = (15, 6)
    _, (ax1, ax3, ax4) = plt.subplots(1, 3)  # Creating plot figure and axes
    ax2 = ax1.twinx()  # Creating separate axis

    # Plotting our LMPs and the given data
    ax1.plot(tP, P/1000, 'k-', label='Pressure Best Fit')
    ax1.plot(tPe, Pe/1000, 'k.', label='Data')
    ax2.plot(tT, T-273.15, 'r-', label='Temperature Best Fit')
    ax2.plot(tTe, Te-273.15, 'r.', label='Data')

    # Plotting misfits
    lim = np.full(len(tT),0)
    ax3.plot(tT, lim, 'k--')
    ax4.plot(tP, lim, 'k--')
    ax3.plot(tT, Tmf, 'rx', label='Temperature Misfit')
    ax4.plot(tP, Pmf/1000, 'kx', label='Pressure Misfit')

    # Setting y limits for each axes, drawing labels and legends
    ax1.set_ylabel('Pressure (kPa)')
    ax2.set_ylabel('Temperature ($^{0}C$)')
    ax3.set_ylabel('Temperature Misfit ($^{0}C$)')
    ax4.set_ylabel('Pressure Misfit (kPa)')
    ax1.set_xlabel('Time (days)')
    ax1.set_title('Pressure and Temperature Models')
    ax1.legend(loc='upper right')
    ax2.legend(loc='lower left')
    ax3.set_title('Temperature Misfit')
    ax4.set_title('Pressure Misfit')

    plt.tight_layout()

    # Either show the plot to the screen or save a version of it to the disk
    os.chdir("../plots")
    save_figure = False

    if not save_figure:
        plt.show()
    else:
        plt.savefig('Pressure_Temp.png', dpi=300)

    return


def temp_forecast():
    ''' Plot our model prediction for temperature over the next 300 days using different injection rates

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":

        It contains commands to plot predictions of temerature change over time using four different injection rates,
        and then display the plot to the screen.

    '''

    # Creating data containing variables necessary for temperature prediction
    data = formulate_models()
    t, P, T, te, Te  = data[0], data[3], data[6], data[7], data[8], 
    Pa = 1291760
    Ta = 453.848

    # Setting initial temperature and pressure values
    T0 = T[-1]
    P0 = P[-1]

    # Allocating time arrays
    t0 = 221
    t1 = 521
    dt = 1
    iterations = int(np.ceil((t1-t0)/dt))
    tp = t0 + np.arange(iterations+1)*dt

    # Creating q1 arrays to interpolate steam injection for 2 full cycles at different rates (0 tonnes/day, 250 tonnes/day, 460 tonnes/day, 1000 tonnes/day)
    q1_0 = np.full(300,0)

    q1_250 = np.full(60,250000)
    q1_250 = np.append(q1_250,np.full(90,0))
    q1_250 = np.append(q1_250,np.full(60,250000))
    q1_250 = np.append(q1_250,np.full(90,0))
   
    q1_460 = np.full(60,460000)
    q1_460 = np.append(q1_460,np.full(90,0))
    q1_460 = np.append(q1_460,np.full(60,460000))
    q1_460 = np.append(q1_460,np.full(90,0))
    
    q1_1000 = np.full(60,1000000)
    q1_1000 = np.append(q1_1000,np.full(90,0))
    q1_1000 = np.append(q1_1000,np.full(60,1000000))
    q1_1000 = np.append(q1_1000,np.full(90,0))
    
    # Creating q2 arrays to interpolate oil/water extraction for 2 full cycles at different injection rates (0 tonnes/day, 250 tonnes/day, 460 tonnes/day, 1000 tonnes/day)
    q2 = interpolate_mass_sink(t)

    q2_0 = np.full(300,0)

    q2_250 = np.full(60,0)
    q2_250 = np.append(q2_250,np.full(90,(max(q2)/2)))
    q2_250 = np.append(q2_250,np.full(60,0))
    q2_250 = np.append(q2_250,np.full(90,(max(q2)/2)))

    q2_460 = np.full(60,0)
    q2_460 = np.append(q2_460,np.full(90,max(q2)))
    q2_460 = np.append(q2_460,np.full(60,0))
    q2_460 = np.append(q2_460,np.full(90,max(q2)))

    q2_1000 = np.full(60,0)
    q2_1000 = np.append(q2_1000,np.full(90,max(q2)*2))
    q2_1000 = np.append(q2_1000,np.full(60,0))
    q2_1000 = np.append(q2_1000,np.full(90,max(q2)*2))

    # Setting fitted parameters from past data
    a = 0.16229278 
    b = 0.02694718
    bt = 8.54611944e-01 
    M0 = 4.45931451e+06

    # Initialising parameter arrays
    pars_P = [a, b, Pa]
    pars_T = [a, b, bt, Pa, Ta, M0]

    # Forecasts for different levels of steam injection
    _, Pzero = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_0, q2_0, P0, pars_P)
    _, P250 = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_250, q2_250, P0, pars_P)
    _, P500 = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_460, q2_460, P0, pars_P)
    _, P1000 = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_1000, q2_1000, P0, pars_P)

    tzero, Tzero = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_0, T0, Pzero, pars_T)
    t250, T250 = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_250, T0, P250, pars_T)
    t460, T460 = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_460, T0, P500, pars_T)
    t1000, T1000 = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_1000, T0, P1000, pars_T)

    # Limit of 240 degrees for dissociation of contaminants
    tlim = np.concatenate((t,tp))
    lim = np.full(len(tlim),240)

    # Creating plot figure and axes
    plt.rcParams["figure.figsize"] = (11, 7)
    f, ax1 = plt.subplots(1, 1) 

    # Plotting our LPMs and the given data
    a, = ax1.plot(t, T-273.15, 'k-', label='Best Fit')
    b, = ax1.plot(te, Te-273.15, 'k.', label='Data')
    c, = ax1.plot(tlim, lim, 'k--', label='Dissociation of Toxic Contaminants')

    d, = ax1.plot(tzero, Tzero-273.15, 'crimson', label='0 tonnes/day')
    e, = ax1.plot(t250, T250-273.15, 'springgreen', label='250 tonnes/day')
    f, = ax1.plot(t460, T460-273.15, 'royalblue', label='460 tonnes/day')
    g, = ax1.plot(t1000, T1000-273.15, 'darkviolet', label='1000 tonnes/day')

    # Drawing labels and legends
    ax1.set_ylabel('Temperature ($^{0}C$)')
    ax1.set_xlabel('Time (days)')
    ax1.set_title('Temperature Forecast')

    # Creating two separate legends for past and forecast data
    first_legend = plt.legend(handles=[a,b], title="Historical Model", loc='upper left')
    second_legend = plt.legend(handles=[c], loc='upper right')
    plt.legend(handles=[d,e,f,g], title="Forecast Models", loc='lower left')
    plt.gca().add_artist(first_legend)
    plt.gca().add_artist(second_legend)

    plt.show()

    return

def uncertainty():
    ''' Plot the final image of Temerature over time using four different injection rates, with uncertainty.

        Parameters:
        -----------
        none

        Returns:
        --------
        none

        Notes:
        ------
        This function called within if __name__ == "__main__":

    '''

    # Creating data containing variables necessary for temperature prediction
    data = formulate_models()
    t, P, T, te, Te  = data[0], data[3], data[6], data[7], data[8], 
    Pa = 1291760
    Ta = 453.848

    # Setting initial temperature and pressure values
    T0 = T[-1]
    P0 = P[-1]

    # Allocating time arrays
    t0 = 221
    t1 = 521
    dt = 1
    iterations = int(np.ceil((t1-t0)/dt))
    tp = t0 + np.arange(iterations+1)*dt

    # Creating q1 arrays to interpolate steam injection for 2 full cycles at different rates (0 tonnes/day, 250 tonnes/day, 460 tonnes/day, 1000 tonnes/day)
    q1_0 = np.full(300,0)

    q1_250 = np.full(60,250000)
    q1_250 = np.append(q1_250,np.full(90,0))
    q1_250 = np.append(q1_250,np.full(60,250000))
    q1_250 = np.append(q1_250,np.full(90,0))
   
    q1_460 = np.full(60,460000)
    q1_460 = np.append(q1_460,np.full(90,0))
    q1_460 = np.append(q1_460,np.full(60,460000))
    q1_460 = np.append(q1_460,np.full(90,0))
    
    q1_1000 = np.full(60,1000000)
    q1_1000 = np.append(q1_1000,np.full(90,0))
    q1_1000 = np.append(q1_1000,np.full(60,1000000))
    q1_1000 = np.append(q1_1000,np.full(90,0))
    
    # Creating q2 arrays to interpolate oil/water extraction for 2 full cycles at different injection rates (0 tonnes/day, 250 tonnes/day, 460 tonnes/day, 1000 tonnes/day)
    q2 = interpolate_mass_sink(t)

    q2_0 = np.full(300,0)

    q2_250 = np.full(60,0)
    q2_250 = np.append(q2_250,np.full(90,(max(q2)/2)))
    q2_250 = np.append(q2_250,np.full(60,0))
    q2_250 = np.append(q2_250,np.full(90,(max(q2)/2)))

    q2_460 = np.full(60,0)
    q2_460 = np.append(q2_460,np.full(90,max(q2)))
    q2_460 = np.append(q2_460,np.full(60,0))
    q2_460 = np.append(q2_460,np.full(90,max(q2)))

    q2_1000 = np.full(60,0)
    q2_1000 = np.append(q2_1000,np.full(90,max(q2)*2))
    q2_1000 = np.append(q2_1000,np.full(60,0))
    q2_1000 = np.append(q2_1000,np.full(90,max(q2)*2))

    # Parameter values
    a = 0.16229278 
    b = 0.02694718
    bt = 8.54611944e-01 
    M0 = 4.45931451e+06

    # Calculating variance for each parameter
    var = 10
    a_var = a/var
    b_var = b/var
    bt_var = bt/var
    M0_var = M0/var

    # Generating normal distributions for parameters
    a_norm = np.random.normal(a,a_var,1000)
    b_norm = np.random.normal(b,b_var,1000)
    bt_norm = np.random.normal(bt,bt_var,1000)
    M0_norm = np.random.normal(M0,M0_var,1000)

    # Creating plot figure and axes
    plt.rcParams["figure.figsize"] = (11, 7)
    f, ax1 = plt.subplots(1, 1) 

    for i in range(200):

        pars_P = [a_norm[i], b_norm[i], Pa]
        pars_T = [a_norm[i], b_norm[i], bt_norm[i], Pa, Ta, M0_norm[i]]

        # Forecasts for different levels of steam injection
        _, Pzero = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_0, q2_0, P0, pars_P)
        _, P250 = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_250, q2_250, P0, pars_P)
        _, P500 = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_460, q2_460, P0, pars_P)
        _, P1000 = solve_ode_pressure(pressure_ode_model, t0, t1, dt, q1_1000, q2_1000, P0, pars_P)

        tzero, Tzero = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_0, T0, Pzero, pars_T)
        t250, T250 = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_250, T0, P250, pars_T)
        t460, T460 = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_460, T0, P500, pars_T)
        t1000, T1000 = solve_ode_temp(temp_ode_model, t0, t1, dt, q1_1000, T0, P1000, pars_T)

        # Limit of 240 degrees for dissociation of contaminants
        tlim = np.concatenate((t,tp))
        lim = np.full(len(tlim),240)

        # Plotting our LPMs and the given data
        a, = ax1.plot(t, T-273.15, 'k-', label='Best Fit')
        b, = ax1.plot(te, Te-273.15, 'k.', label='Data')
        c, = ax1.plot(tlim, lim, 'k--', label='Dissociation of Toxic Contaminants')

        d, = ax1.plot(tzero, Tzero-273.15, color="crimson", label='0 tonnes/day',linewidth=0.3, alpha=0.7)
        e, = ax1.plot(t250, T250-273.15 ,color="springgreen", label='250 tonnes/day',linewidth=0.3, alpha=0.7)
        f, = ax1.plot(t460, T460-273.15,color="royalblue", label='460 tonnes/day',linewidth=0.3, alpha=0.7)
        g, = ax1.plot(t1000, T1000-273.15,color="darkviolet", label='1000 tonnes/day',linewidth=0.3, alpha=0.7)

        i = i + 1

    # Drawing labels and legends
    ax1.set_ylabel('Temperature ($^{0}C$)')
    ax1.set_xlabel('Time (days)')
    ax1.set_title('Temperature Forecast')

    # Creating two separate legends for past and forecast data
    first_legend = plt.legend(handles=[a,b], title="Historical Model", loc='upper left')
    second_legend = plt.legend(handles=[c], loc='upper right')
    plt.legend(handles=[d,e,f,g], title="Forecast Models", loc='lower left')
    plt.gca().add_artist(first_legend)
    plt.gca().add_artist(second_legend)

    plt.show()

if __name__ == "__main__":
    #plot_models()
    #temp_forecast()
    uncertainty()

    
