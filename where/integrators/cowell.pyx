"""Integration by Cowell method.

Description:
------------

Integrating equation  d/dt (state) = integrand(state, t)

with initial value integrand(0) = initial_state


References:

[1] Oesterwinter, C. & Cohen, C.J. Celestial Mechanics (1972) 5: 317. https://doi.org/10.1007/BF01228434
    page 329 - 333
[2] Montenbruck, Oliver and Gill, Eberhard: Satellite Orbits, Springer Verlag, 2000.




$Revision$
$Date$
$LastChangedBy$
"""

#Midgard imports
from midgard.dev import log

#Standard imports
import numpy as np
import math
import scipy.special

def integrate(integrand, double[:] initial_state, double grid_step, double end_time):
    """
    Input:
        integrand:        Function f in the equation dy/dt = f(y,t) to be solved
        initial_state:    y at time t = 0.
        grid_step:        step length in integrator
        end_time:         End time of integration

    Output:
        orbit:          Table of values of y at times in time_grid
        time_grid:      Times corresponding to y above
    """
    # For consistency with the notation in the equations in [1], we set
    cdef double h = grid_step
    cdef int idx, i, j
    cdef int n = 0 
    # Coefficients from page 329 - 332 in [1].
    cdef int a = 7
    cdef int b = 7
    cdef int c = 14
    cdef int m = 12
    cdef double D = 32011868528640000

    cdef long[:, :] ac_vector = get_ac()
    cdef long[:] a_vector = ac_vector[0, :] 
    cdef long[:] c_vector = ac_vector[1, :]
    cdef double[:] gamma_prime = compute_gamma_prime(a_vector, c_vector, a, b)
    cdef double[:] gamma = compute_gamma(gamma_prime, a, b)
    cdef double[:] b_vector = compute_b(a, b)
    cdef double[:] alpha_star = compute_alpha_star(c, c_vector)
    cdef double[:] beta_star = compute_beta_star(c, a_vector)
    cdef double[:] c_prime = compute_c_prime(m, c_vector)
    cdef double[:] alpha_prime = compute_alpha_prime(m, c_prime)
    cdef double[:] a_prime = compute_a_prime(m, a_vector)
    cdef double[:] beta_prime = compute_beta_prime(m, a_prime)
    cdef double[:] alpha = compute_alpha(m, c_vector)
    cdef double[:] beta = compute_beta(m, a_vector)
    
    # Set up time grid
    cdef double[:] time_grid = np.arange(-c * h, math.ceil(end_time / h) * h + h, h)
    cdef int num_steps = math.ceil(end_time / h) + c + 1

    cdef double[:, :] r = np.zeros((num_steps, 3))
    cdef double[:, :] v = np.zeros((num_steps, 3))
    cdef double[:, :] acc = np.zeros((num_steps, 3))

    # The 6x6 transition matrix phi consists of phi_r (upper 3x6 part) and phi_v (lower 3x6 part),
    # see [2], section 7.2
    cdef double[:, :] phi_r = np.zeros((num_steps, 18))
    cdef double[:, :] phi_v = np.zeros((num_steps, 18))
    cdef double[:, :] phi_acc = np.zeros((num_steps, 18))

    # The 6xnum_params matrix S consists of s_r (upper 3xnum_params part) and s_v (lower 3xnum_params part)
    cdef int num_params = int((initial_state.shape[0] - 42) / 6)

    cdef double[:, :] s_r = np.zeros((num_steps, 3 * num_params))
    cdef double[:, :] s_v = np.zeros((num_steps, 3 * num_params))
    cdef double[:, :] s_acc = np.zeros((num_steps, 3 * num_params))

    for idx in range(0, 3):
        r[0, idx] = initial_state[idx]
        v[0, idx] = initial_state[3 + idx]
        acc[0, idx] = integrand(initial_state, n)[3 + idx]
    for idx in range(0, 18):
        phi_r[0, idx] = initial_state[6 + idx]
        phi_v[0, idx] = initial_state[24 + idx]
        phi_acc[0, idx] = integrand(initial_state, n)[24 + idx]
    for idx in range(0, 3 * num_params):
        s_r[0, idx] = initial_state[42 + idx]
        s_v[0, idx] = initial_state[42 + 3 * num_params + idx]
        s_acc[0, idx] = integrand(initial_state, n)[42 + 3 * num_params + idx]

    # Set accelerations (and partials) equal to initial value in the range from -a to b
    for n in range(-a, b + 1):
        for idx in range(0, 3):
            acc[n, idx] = acc[0, idx]
        for idx in range(0, 18):
            phi_acc[n, idx] = phi_acc[0, idx]
        for idx in range(0, 3 * num_params):
            s_acc[n, idx] = s_acc[0, idx]

    # Setting up starting table
    # TODO: Add convergence check (p 332)
    for K in range(0, 4):
        # Equation (6) in [1]
        for idx in range(0, 3):
            r[-1, idx] = r[0, idx] - h * v[0, idx]
        for idx in range(0, 18):
            phi_r[-1, idx] = phi_r[0, idx] - h * phi_v[0, idx]
        for idx in range(0, 3 * num_params):
            s_r[-1, idx] = s_r[0, idx] - h * s_v[0, idx]
        for i in range(0, a + b + 1):
            for idx in range(0, 3):
                r[-1, idx] += h**2 * gamma[i] * acc[b - i, idx] / D
            for idx in range(0, 18):
                phi_r[-1, idx] += h**2 * gamma[i] * phi_acc[b - i, idx] / D
            for idx in range(0, 3 * num_params):
                s_r[-1, idx] += h**2 * gamma[i] * s_acc[b - i, idx] /D

        # Extrapolation beyond -a and b to c lines before and after epoch by
        for n in range(a + 1, c + 1):
            for idx in range(0, 3):
                acc[-n, idx] = 0
            for idx in range(0, 18):
                phi_acc[ -n, idx] = 0
            for idx in range(0, 3 * num_params):
                s_acc[-n, idx] = 0

            for j in range(0, a + b + 1):
                for idx in range(0, 3):
                    acc[-n, idx] += b_vector[j] * acc[-n + 1 + j, idx]
                for idx in range(0, 18):
                    phi_acc[-n, idx] += b_vector[j] * phi_acc[-n + 1 + j, idx]
                for idx in range(0, 3 * num_params):
                    s_acc[-n, idx] += b_vector[j] * s_acc[-n + 1 + j, idx]

        for n in range(b + 1, c + 1):
            for idx in range(0, 3):
                acc[n, idx] = 0
            for idx in range(0, 18):
                phi_acc[n, idx] = 0
            for idx in range(0, 3 * num_params):
                s_acc[n, idx] = 0

            for j in range(0, a + b + 1):
                for idx in range(0, 3):
                    acc[n, idx] += b_vector[j] * acc[n - 1 - j, idx]
                for idx in range(0, 18):
                    phi_acc[n, idx] += b_vector[j] * phi_acc[n - 1 - j, idx]
                for idx in range(0, 3 * num_params):
                    s_acc[n, idx] += b_vector[j] * s_acc[n - 1 - j, idx]

        # Position and velocity in the range -a to b
        for n in range(1, b + 1):
            for idx in range(0, 3):
                r[n, idx] = 2 * r[n - 1, idx] - r[n - 2, idx]
                v[n, idx] = v[n - 1, idx]
            for idx in range(0, 18):
                phi_r[n, idx] = 2 * phi_r[n - 1, idx] - phi_r[n - 2, idx]
                phi_v[n, idx] = phi_v[n - 1, idx]
            for idx in range(0, 3 * num_params):
                s_r[n, idx] = 2 * s_r[n - 1, idx] - s_r[n - 2, idx]
                s_v[n, idx] = s_v[n - 1, idx]

            for j in range(0, c + 1):
                for idx in range(0, 3):
                    r[n, idx] += h**2 * alpha_star[j] * acc[n - j, idx] / D
                    v[n, idx] += h * beta_star[j] * acc[n - j, idx] / D
                for idx in range(0, 18):
                    phi_r[n, idx] += h**2 * alpha_star[j] * phi_acc[n - j, idx] / D
                    phi_v[n, idx] += h * beta_star[j] * phi_acc[n - j, idx] / D
                for idx in range(0, 3 * num_params):
                    s_r[n, idx] += h**2 * alpha_star[j] * s_acc[n - j, idx] / D
                    s_v[n, idx] += h * beta_star[j] * s_acc[n - j, idx] / D

        for n in range(2, a + 1):
            for idx in range(0, 3):
                r[-n, idx] = 2 * r[-n + 1, idx] - r[-n + 2, idx]
            for idx in range(0, 18):
                phi_r[-n, idx] = 2 * phi_r[-n + 1, idx] - phi_r[-n + 2, idx]
            for idx in range(0, 3 * num_params):
                s_r[-n, idx] = 2 * s_r[-n + 1, idx] - s_r[-n + 2, idx]

            for j in range(0, c + 1):
                for idx in range(0, 3):
                    r[-n, idx] += h**2 * alpha_star[j] * acc[-n + j, idx] / D
                for idx in range(0, 18):
                    phi_r[-n, idx] += h**2 * alpha_star[j] * phi_acc[-n + j, idx] / D
                for idx in range(0, 3 * num_params):
                    s_r[-n, idx] += h**2 * alpha_star[j] * s_acc[-n + j, idx] / D

        for n in range(1, a + 1):
            for idx in range(0, 3):
                v[-n, idx] = v[-n + 1, idx]
            for idx in range(0, 18):
                phi_v[-n, idx] = phi_v[-n + 1, idx]
            for idx in range(0, 3 * num_params):
                s_v[-n, idx] = s_v[-n + 1, idx]

            for j in range(0, c + 1):
                for idx in range(0, 3):
                    v[-n, idx] += -h * beta_star[j] * acc[-n + j, idx] / D
                for idx in range(0, 18):
                    phi_v[-n, idx] += -h * beta_star[j] * phi_acc[-n + j, idx] / D
                for idx in range(0, 3 * num_params):
                    s_v[-n, idx] += -h * beta_star[j] * s_acc[-n + j, idx] / D

        for n in range(-a, b + 1):
            state_vector = integrand(np.hstack((r[n], v[n], phi_r[n], phi_v[n], s_r[n], s_v[n])), n)
            for idx in range(0, 3):
                acc[n, idx] = state_vector[3 + idx]
            for idx in range(0, 18):
                phi_acc[n, idx] = state_vector[24 + idx]
            for idx in range(0, 3 * num_params):
                s_acc[n, idx] = state_vector[42 + 3 * num_params + idx]


    # Next step after starting routine
    for n in range(b + 1, math.ceil(end_time/h) + 1):
        for idx in range(0, 3):
            r[n, idx] = 2 * r[n - 1, idx] - r[n - 2, idx]
            v[n, idx] = v[n - 1, idx]
        for idx in range(0, 18):
            phi_r[n, idx] = 2 * phi_r[n - 1, idx] - phi_r[n - 2, idx]
            phi_v[n, idx] = phi_v[n - 1, idx]
        for idx in range(0, 3 * num_params):
            s_r[n, idx] = 2 * s_r[n - 1, idx] - s_r[n - 2, idx]
            s_v[n, idx] = s_v[n - 1, idx]

        for j in range(0, m + 1):
            for idx in range(0, 3):
                r[n, idx] += h**2 * alpha_prime[j] * acc[n - 1 - j, idx] / D
                v[n, idx] += h * beta_prime[j] * acc[n - 1 - j, idx] / D
            for idx in range(0, 18):
                phi_r[n, idx] += h**2 * alpha_prime[j] * phi_acc[n - 1 - j, idx] / D
                phi_v[n, idx] += h * beta_prime[j] * phi_acc[n - 1 - j, idx] / D
            for idx in range(0, 3 * num_params):
                s_r[n, idx] += h**2 * alpha_prime[j] * s_acc[n - 1 - j, idx] / D
                s_v[n, idx] += h * beta_prime[j] * s_acc[n - 1 - j, idx] /D

        state_vector = integrand(np.hstack((r[n], v[n], phi_r[n], phi_v[n], s_r[n], s_v[n])), n)
        for idx in range(0, 3):
            acc[n, idx] = state_vector[3 + idx]
        for idx in range(0, 18):
            phi_acc[n, idx] = state_vector[24 + idx]
        for idx in range(0, 3 * num_params):
            s_acc[n, idx] = state_vector[42 + num_params * 3 + idx]

        # Coordinates are refined with corrector formulae
        for idx in range(0, 3):
            r[n, idx] = 2 * r[n - 1, idx] - r[n - 2, idx]
            v[n, idx] = v[n - 1, idx]
        for idx in range(0, 18):
            phi_r[n, idx] = 2 * phi_r[n - 1, idx] - phi_r[n - 2, idx]
            phi_v[n, idx] = phi_v[n - 1, idx]
        for idx in range(0, 3 * num_params):
            s_r[n, idx] = 2 * s_r[n - 1, idx] - s_r[n - 2, idx]
            s_v[n, idx] = s_v[n - 1, idx]

        for j in range(0, m + 1):
            for idx in range(0, 3):
                r[n, idx] += h**2 * alpha[j] * acc[n - j, idx] / D
                v[n, idx] += h * beta[j] * acc[n - j, idx] / D
            for idx in range(0, 18):
                phi_r[n, idx] += h**2 * alpha[j] * phi_acc[n - j, idx] / D
                phi_v[n, idx] += h * beta[j] * phi_acc[n - j, idx] / D
            for idx in range(0, 3 * num_params):
                s_r[n, idx] += h**2 * alpha[j] * s_acc[n - j, idx] / D
                s_v[n, idx] += h * beta[j] * s_acc[n - j, idx] / D

        state_vector = integrand(np.hstack((r[n], v[n], phi_r[n], phi_v[n], s_r[n], s_v[n])), n)
        for idx in range(0, 3):
            acc[n, idx] = state_vector[3 + idx]
        for idx in range(0, 18):
            phi_acc[n, idx] = state_vector[24 + idx]
        for idx in range(0, num_params * 3):
            s_acc[n, idx] = state_vector[42 + num_params * 3 + idx]

    orbit = np.hstack((r, v, phi_r, phi_v, s_r, s_v))

    return orbit[:-c, :], time_grid[c:]

cdef long[:, :] get_ac():
    """Coefficients from page 329 in [1]
    """
    cdef long[:, :] ac
    ac = np.array([[32011868528640000, -16005934264320000, -2667655710720000, -1333827855360000, -844757641728000,
                  -600222534912000, -456783110784000, -363891528000000, -299520219398400, -252655401398400,
                  -217227737563200, -189640115028000, -167636336098320, -149735464049160, -134928496929540,
                  -122506205369730, -111956703448001], [32011868528640000, -32011868528640000, 2667655710720000, 0, 
                  -133382785536000, -133382785536000, -116974585728000, -100566385920000, -86707632211200, 
                  -75398324601600, -66193573118400, -58648487788800, -52401453198480, -47174128491600, -42755108505900, 
                  -38983584907800, -35736323456205]])

    return ac


cdef double[:] compute_gamma_prime(a_vector, c_vector, a, b):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] gamma_prime = np.zeros(a + b + 1)
    for j in range(0, a + b + 1):
        for i in range(0, min(j + 1, b) + 1):
            gamma_prime[j] += (-1)**i * scipy.special.binom(b, i) * (a_vector[j + 1 - i] - c_vector[j + 1 - i])
    return gamma_prime

cdef double[:] compute_gamma(gamma_prime, a, b):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] gamma = np.zeros(a + b + 1)
    for i in range(0, a + b + 1):
        for j in range(i, a + b + 1):
            gamma[i] += scipy.special.binom(j, i) * gamma_prime[j]
        gamma[i] *= (-1)**i
    return gamma

cdef double[:] compute_b(a, b):
    """Method from page 330 in [1]
    """
    cdef int j
    cdef double[:] b_vector = np.zeros(a + b + 1)
    for j in range(0, a + b + 1):
        b_vector[j] = (-1)**j * scipy.special.binom(a + b + 1, j + 1)
    return b_vector

cdef double[:] compute_alpha_star(c, c_vector):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] alpha_star = np.zeros(c + 1)
    for j in range(0, c + 1):
        for i in range(j, c + 1):
            alpha_star[j] += (-1)**j * scipy.special.binom(i, j) * c_vector[i]
    return alpha_star

cdef double[:] compute_beta_star(c, a_vector):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] beta_star = np.zeros(c + 1)
    for j in range(0, c + 1):
        for i in range(j, c + 1):
            beta_star[j] += (-1)**j * scipy.special.binom(i, j) * a_vector[i]
    return beta_star

cdef double[:] compute_c_prime(m, c_vector):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] c_prime = np.zeros(m + 1)
    for i in range(0, m + 1):
        for j in range(0, i + 1):
            c_prime[i] += c_vector[j]
    return c_prime

cdef double[:] compute_alpha_prime(m, c_prime):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] alpha_prime = np.zeros(m + 1)
    for j in range(0, m + 1):
        for i in range(j, m + 1):
            alpha_prime[j] += (-1)**j * scipy.special.binom(i, j) * c_prime[i]
    return alpha_prime

cdef double[:] compute_a_prime(m, a_vector):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] a_prime = np.zeros(m + 1)
    for i in range(0, m + 1):
        for j in range(0, i + 1):
            a_prime[i] += a_vector[j]
    return a_prime

cdef double[:] compute_beta_prime(m, a_prime):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] beta_prime = np.zeros(m + 1)
    for j in range(0, m + 1):
        for i in range(j, m + 1):
            beta_prime[j] += (-1)**j * scipy.special.binom(i, j) * a_prime[i]
    return beta_prime

cdef double[:] compute_alpha(m, c_vector):
    """Method from page 330 in [1]
    """
    cdef int i, j
    cdef double[:] alpha = np.zeros(m + 1)
    for j in range(0, m + 1):
        for i in range(j, m + 1):
            alpha[j] += (-1)**j * scipy.special.binom(i, j) * c_vector[i]
    return alpha

cdef double[:] compute_beta(m, a_vector):
    """Method from page 330 in [1]
    """
    cdef int i, j
    beta = np.zeros(m + 1)
    for j in range(0, m + 1):
        for i in range(j, m + 1):
            beta[j] += (-1)**j * scipy.special.binom(i, j) * a_vector[i]
    return beta
