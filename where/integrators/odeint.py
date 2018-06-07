"""Integration by scipy package odeint

Description:
------------

Integrating equation  d/dt (state) = integrand(state, t)

with initial value integrand(0) = initial_state




"""

# Where imports
from where.lib import plugins

# Standard imports
import numpy as np
from scipy import integrate


@plugins.register
def odeint(integrand, initial_state, grid_step, end_time):

    # Set up time grid
    time_grid = np.arange(0, end_time + grid_step, grid_step)

    # Integrate orbit
    orbit = integrate.odeint(integrand, initial_state, time_grid, hmax=grid_step)

    return orbit, time_grid
