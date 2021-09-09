import numpy as np
from main import solve_ode_pressure, solve_ode_temp, temp_ode_model, pressure_ode_model
from numpy.linalg import norm


def test_solve_ode_pressure(f, t0, t1, dt, q1, q2, P0, pars, t_soln, P_soln):
   
    '''
    This function tests the Improved Euler solver function for pressure.
    '''
    t, P = solve_ode_pressure(f, t0, t1, dt, q1, q2, P0, pars)
    
    # Checks that the returned output from the solver function and the solution 
    # are the same, with an error margin of 1e-10
    assert norm(t - t_soln) < 1e-10
    assert norm(P - P_soln) < 1e-10


def test_solve_ode_temp(f, t0, t1, dt, q1, T0, P, pars, t_soln, T_soln):

    '''
    This function tests the Improved Euler solver function for temperature.
    '''
    t, T = solve_ode_temp(f, t0, t1, dt, q1, T0, P, pars)

    # Checks that the returned output from the solver function and the solution 
    # are the same, with an error margin of 1e-10
    assert norm(t - t_soln) < 1e-10
    assert norm(T - T_soln) < 1e-10

def test_solve_ode_pressure_0(f, t0, t1, dt, q1, q2, P0, pars, t_soln, P_soln):
   
    '''
    This function tests the Improved Euler solver function for pressure.
    It tests that a ZeroDivisionError is raised when input dt = 0
    '''
    try:
        t, P = solve_ode_pressure(f, t0, t1, dt, q1, q2, P0, pars)
    except ZeroDivisionError:
        assert True
    
def test_solve_ode_temp_0(f, t0, t1, dt, q1, T0, P, pars, t_soln, T_soln):

    '''
    This function tests the Improved Euler solver function for temperature.
    It tests that a ZeroDivisionError is raised when M0, a or dt are zero.
    '''
    try:
        t, T = solve_ode_temp(f, t0, t1, dt, q1, T0, P, pars)
    except ZeroDivisionError:
        assert True


if __name__ == "__main__":
    
    test_solve_ode_pressure(f=pressure_ode_model, t0=0, t1=1, dt=1, q1=np.array([4]), q2=np.array([4]), P0=2, pars=[2,2,2], t_soln=np.array([0,1]), P_soln=np.array([2,4.25]))   
    
    test_solve_ode_pressure_0(f=pressure_ode_model, t0=0, t1=2, dt=0, q1=np.array([4]), q2=np.array([7]), P0=8, pars=[2,2,2], t_soln=0, P_soln=8)
    
    test_solve_ode_temp(f=temp_ode_model, t0=0, t1=1, dt=1, q1=np.array([4]), T0=500, P=np.array([5]), pars=[1, 2, 5, 3, 2, 4], t_soln=np.array([0,1]), T_soln=np.array([500, 1956.85]))

    test_solve_ode_temp_0(f=temp_ode_model, t0=0, t1=2, dt=0, q1=np.array([4]), T0=500, P=np.array([5]), pars=[0, 2, 5, 3, 2, 0], t_soln=0, T_soln=500)