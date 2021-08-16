#imports
import numpy as np
from numpy.linalg import norm, solve
from matplotlib import pyplot as plt

def solve_ode(f, t0, t1, dt, x0, pars):

	''' Solve an ODE numerically.

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
		x0 : float
			Initial value of solution.
		pars : array-like
			List of parameters passed to ODE function f.

		Returns:
		--------
		t : array-like
			Independent variable solution vector.
		x : array-like
			Dependent variable solution vector.

		Notes:
		------
		ODE should be solved using the Improved Euler Method. 

		Function q(t) should be hard coded within this method. Create duplicates of 
		solve_ode for models with different q(t).

		Assume that ODE function f takes the following inputs, in order:
			1. independent variable
			2. dependent variable
			3. forcing term, q
			4. all other parameters
	'''
	
	iterations = int(np.ceil((t1-t0)/dt))		# compute number of Euler steps to take
	t = t0 + np.arange(iterations+1)*dt			# x array
	x = 0.*t							# array to store solution
	x[0] = x0							# set initial value
	q = -1
	
	for i in range(iterations):
		# Creating a loop to find the next euler step and iterate the next solution
		k1 = f(t[i], x[i], q, *pars)
		k2 = f(t[i]+dt, x[i]+dt*k1, q, *pars)
        x[i+1] = x[i] + 0.5*dt*(k1+k2)

	return t, x

def temperature_ode_model(q,T,m,P,a,b,d,T0,P0):
	
	''' Create a model for the bitumen production.
    
    	Parameters:
    	-----------
        t : float
            Independent variable.
        x : float
            Dependent variable.
        q : float
            Source/sink rate.
        a : float
            Source/sink strength parameter.
        b : float
            Recharge strength parameter.
        x0 : float
            Ambient value of dependent variable.
    
		Returns:
    	--------
    	dTdt : float
        The rate of change in temperature with respect to time

    '''
	
	dTdt = q*T/m + (T*b)/(m*a)*(P-P0) + d*(T-T0)
	#returning variable
	return dTdt