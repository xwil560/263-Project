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
	# initialise
	nt = int(np.ceil((t1-t0)/dt))		# compute number of Euler steps to take
	ts = t0+np.arange(nt+1)*dt			# x array
	xs = 0.*ts							# array to store solution
	xs[0] = x0							# set initial value
	for i in range(1,nt+1):
		#create a loop to find the next euler step and iterate the next solution
		dxdtk=f(ts[i-1],xs[i-1],*pars)
		dxdtk2=f(ts[i-1]+dt,xs[i-1]+dxdtk*dt,*pars)
		xs[i]=xs[i-1]+dt/2*(dxdtk+dxdtk2)

	return ts,xs
	#pass
def LMP(q,T,m,p,a,b,d,T0,P0):
    """ Create a model for the bitumen production.
    
    Parameters:
    -----------
    
    Returns:
    --------
    dTdt: float
        the rate of change in temperature
    """
    return q*T/m+(T*b)/(m*a)*(P-P0)+d*(T-T0)
    
if __name__=="__main__":


    