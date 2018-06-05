#!/usr/bin/env python3
"""Where library module additional mathematical functions

Example:
--------

    from where.lib import mathp
    ...

Description:
------------

This module will provide functions additional mathematical functions.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

# External library imports
from scipy import interpolate
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt

# Where imports
from where.lib import log


def _determine_time_derivate(func, x_new, num_col, t_offset=0.5):
    """Determine time derivate based on interpolated values.

    Args:
        func  (function)       - Interpolator function
        x_new (numpy.ndarray)  - 1-dimensional array with x-values (e.g. time entries in [s]), for which interpolation
                                 should be done.
        num_col (int)          - Number of original y-array columns
        t_offset (float)       - Time offset used to determine interpolated x-values before and after given x_new
                                 values in [s]

    Returns:
        numpy.ndarray:         - n-dimensional array with time derivate of y-values (e.g. orbit velocity vectors).
    """
    if isinstance(x_new, float):
        num_row = 1
    else:
        num_row = len(x_new)

    y_new_dot = np.zeros((num_row, num_col))

    y1 = func(x_new - t_offset)
    y2 = func(x_new + t_offset)

    if y1.ndim == 1:
        y1 = y1[np.newaxis]
        y2 = y2[np.newaxis]

    if num_col == 1:
        y_new_dot = (y2 - y1) / (t_offset * 2)

    else:
        for idx in range(0, num_col):
            y_new_dot[:, idx] = (y2[:, idx] - y1[:, idx]) / (t_offset * 2)

    return y_new_dot


def _get_nearest_idx_for_given_data(given_data, data):
    """Get nearest indices for each element of a dataset, which corresponds to a given dataset

    This can for example be used to find the nearest observation time entry in the initial precise orbit data.

    Args:
        given_data (numpy.ndarray): original given dataset (e.g. orbit data)
        data (numpy.ndarray): data set (e.g. observation epochs), for which nearest value should be searched in given
                              data
    Returns:
        tuple: with following entries

    ====================  =================  ==========================================================================
     Entry                 Type               Description
    ====================  =================  ==========================================================================
     nearest_idx           list               Indices, which represents the nearest value entry for a given set of data
     is_identical          list               List, where: True indicates that the given and searched value is identical
                                                and False not
    ====================  =================  ==========================================================================
    """
    nearest_idx = list()
    is_identical = list()

    # Loop over all searching data
    for val in data:

        # Pick all 'given_data' values before 'val'
        given_values_less_val = np.less_equal(given_data, val)

        if not np.any(given_values_less_val):
            log.fatal(
                "Given data range (from {} to {}) corresponds not to the data value {}.",
                np.min(given_data),
                np.max(given_data),
                val,
            )

        # Convert True/False array to indices
        idx = given_values_less_val.nonzero()[0]

        # Find nearest index in given data
        nearest_idx.append(idx[np.argmax(given_data[idx])])

        # Check if an identical value can be found in 'given_data'
        if val in given_data:
            is_identical.append(True)
        else:
            is_identical.append(False)

    return nearest_idx, is_identical


def moving_window_interpolation(x, y, x_new, model="lagrange", window_size=10, plot=False):
    """Interpolation by using a moving window

    The first step is to define the moving window for the given 'x' array by looping over all 'x_new' array elements.
    Hereby the nearest/or equal 'x' value is searched for the given 'x_new' value. At the end the array 'x_idx' is
    created, which includes the indices for the nearest/or equal 'x' value. Afterwards the 'x_idx' array indices are
    used to define the moving window size by adding the 'window_size'. Then the needed values for the interpolation are
    selected for the moving window and the interpolation is done for the 'x_new' values.

    Args:
        x (numpy.ndarray)      - 1-dimensional array with original x-values (e.g. time entries in [s]).
        y (numpy.ndarray)      - n-dimensional array with original y-values (e.g. orbit position vectors).
        x_new (numpy.ndarray)  - 1-dimensional array with x-values (e.g. time entries in [s]), for which interpolation
                                 should be done.
        model (str)            - Interpolation model, which can be:
                                    InterpolatedUnivariateSpline
                                    interp1d
                                    lagrange (limited to around 20 data points)
        window_size (int)      - size of moving window, which defines number of datapoints before and after the
                                 interpolation entry
        plot (bool)            - Flag (True|False) for plotting original data against interpolated ones.

    Returns:
        tuple: with following entries

    ====================  =================  ==========================================================================
     Entry                 Type               Description
    ====================  =================  ==========================================================================
     y_new                 numpy.ndarray      n-dimensional array with interpolated y-values (e.g. orbit position
                                              vectors).
     y_new_derivate        numpy.ndarray      n-dimensional array with time derivate of y-values (e.g. orbit position
                                              vectors).
    ====================  =================  ==========================================================================

    """
    num_x = len(x)
    num_xnew = len(x_new)

    y_new = np.zeros((num_xnew, y.shape[1]))
    y_new_derivate = np.zeros((num_xnew, y.shape[1]))

    # Get belonging indices for 'x_new' values in 'x'
    x_idx, is_identical = _get_nearest_idx_for_given_data(x, x_new)

    # TODO: is_identical array could be used to skip the interpolation, when the interpolated x value is identical with
    #      the given x value. How the y_new_derivate should be determined? scipy.misc.derivate?

    for num_obs, idx in enumerate(x_idx):

        s_idx = idx - window_size + 1  # start index of moving window
        e_idx = idx + window_size + 1  # end index of moving window

        # TODO: would it be better to take at least complete window size amount of data
        if s_idx < 0:
            s_idx = 0  # Index setting at the beginnging of the data series
        if e_idx > num_x:
            e_idx = num_x  # Index setting at the end of the data series

        y_new[num_obs], y_new_derivate[num_obs] = interpolation(
            x[s_idx:e_idx], y[s_idx:e_idx], x_new[num_obs:num_obs + 1], model=model, plot=False
        )

    return y_new, y_new_derivate


def interpolation(x, y, x_new, model="InterpolatedUnivariateSpline", plot=False, title=""):
    """Determination of interpolated data and its derivation by choosing a interpolation method

    Args:
        x (numpy.ndarray)      - 1-dimensional array with original x-values (e.g. time entries in [s]).
        y (numpy.ndarray)      - 2-dimensional array with original y-values (e.g. orbit position vectors).
        x_new (numpy.ndarray)  - 1-dimensional array with x-values (e.g. time entries in [s]), for which interpolation
                                 should be done.
        model (str)            - Interpolation model, which can be:t
                                    InterpolatedUnivariateSpline
                                    interp1d
                                    lagrange (limited to around 20 data points)
        window_size (int)      - size of moving window, which defines number of datapoints before and after the
                                 interpolation entry
        plot (bool)            - Flag (True|False) for plotting original data against interpolated ones.
        title (str)            - Title of plot.

    Returns:
        tuple: with following entries

    ====================  =================  ==========================================================================
     Entry                 Type               Description
    ====================  =================  ==========================================================================
     y_new                 numpy.ndarray      2-dimensional array with interpolated y-values (e.g. orbit position
                                              vectors).
     y_new_derivate        numpy.ndarray      2-dimensional array with time derivate of y-values (e.g. orbit position
                                              vectors).
    ====================  =================  ==========================================================================

    """

    num_col = y.shape[1]

    if isinstance(x_new, float):
        num_row = 1
    else:
        num_row = len(x_new)

    y_new = np.zeros((num_row, num_col))
    y_new_dot = np.zeros((num_row, num_col))

    if model == "InterpolatedUnivariateSpline":

        for idx in range(0, num_col):

            # TODO: InterpolatedUnivariateSpline seems to have problems with multidimensional arrays
            interpolate_ = interpolate.InterpolatedUnivariateSpline(x, y[:, idx], k=3)
            y_new[:, idx] = interpolate_(x_new)
            # y_new_dot[:,idx] = spline.derivative()(x_new) #TODO: Does this work too?
            y_new_dot[:, idx] = _determine_time_derivate(interpolate_, x_new, 1)

        if plot == True:
            _plot_interpolation(x, y, x_new, y_new, title=title)

    elif model == "interp1d":

        # TODO: Extrapolation has to be handled (e.g. with fill_value=(0.0,0.0) argument)
        interpolate_ = interpolate.interp1d(x, y, kind="cubic", axis=0, bounds_error=False)
        y_new = interpolate_(x_new)
        y_new_dot = _determine_time_derivate(interpolate_, x_new, num_col)

        if plot == True:
            _plot_interpolation(x, y, x_new, y_new, title=title)

    elif model == "lagrange":

        for idx in range(0, num_col):

            # Rescaling of data necessary, because Lagrange interpolation is numerical unstable
            xm = np.mean(x)
            ym = np.mean(y[:, idx])
            xs = np.std(x)
            ys = np.std(y[:, idx])
            xscaled = (x - xm) / xs
            yscaled = (y[:, idx] - ym) / ys

            # interpolate_ = interpolate.lagrange(xscaled, yscaled) # worst performance
            # interpolate_ = _lagrange2(xscaled, yscaled) # improved performance
            interpolate_ = _lagrange(xscaled, yscaled)  # fastest performance
            y_new[:, idx] = interpolate_((x_new - xm) / xs) * ys + ym

            # Determine derivate of 'y_new'
            t_offset = 0.5
            y1 = interpolate_((x_new - t_offset - xm) / xs) * ys + ym
            y2 = interpolate_((x_new + t_offset - xm) / xs) * ys + ym

            y_new_dot[:, idx] = (y2 - y1) / (t_offset * 2)

        if plot == True:
            _plot_interpolation(x, y, x_new, y_new, title=title)

    elif model == "BarycentricInterpolator":

        for idx in range(0, num_col):

            # TODO: Is the rescaling of data necessary here?
            xm = np.mean(x)
            ym = np.mean(y[:, idx])
            xs = np.std(x)
            ys = np.std(y[:, idx])
            xscaled = (x - xm) / xs
            yscaled = (y[:, idx] - ym) / ys

            interpolate_ = interpolate.BarycentricInterpolator(xscaled, yscaled)
            y_new[:, idx] = interpolate_((x_new - xm) / xs) * ys + ym

            # Determine derivate of 'y_new'
            t_offset = 0.5
            y1 = interpolate_((x_new - t_offset - xm) / xs) * ys + ym
            y2 = interpolate_((x_new + t_offset - xm) / xs) * ys + ym

            y_new_dot[:, idx] = (y2 - y1) / (t_offset * 2)

        if plot == True:
            _plot_interpolation(x, y, x_new, y_new, title=title)

    #    elif model == 'polyfit':

    #        #y_new_dot = np.zeros((len(x_new), num_col))
    #        degree = 2
    #        poly_dot_coeff = np.zeros((degree, num_col))
    #

    #        #TODO: Polynomial fit over the complete time period of 3 days are not convenient.
    #        #      Difference of up to 1000 km between original SP3 orbit data and polynomial
    #        #      fit.

    #        x_ref = np.amin(x)     #TODO: orb_ref_time should be part of 'orb' Dataset
    #        x = (x - x_ref) * 86400
    #        poly_coeff = poly.polyfit(x, y, degree)
    #        y_new = poly.polyval(x_new,poly_coeff)

    #        for idx in range(0,num_col):
    #            poly_dot_coeff[:,idx] = np.polyder(poly_coeff[:,idx])
    #            polynomial_dot = np.poly1d(poly_dot_coeff[:,idx])
    #            y_new_dot[:,idx] = polynomial_dot((x_new - x_ref) * 86400)

    #        if plot == True:
    #            _plot_orbit_polynomial(sat, poly_coeff, x, y)

    return y_new, y_new_dot


def unit_vector(vector):
    """Determine unit vector based on 1- or 2-dimensional array.

    Args:
        vector (numpy.ndarray)  - 1- or 2-dimensional array with shape (n, m)

    Returns:
        numpy.ndarray: 1-dimensional unit vector of length n
    """
    unit_vector = np.zeros((len(vector), vector.shape[1]))
    norm = np.linalg.norm(vector, axis=1)
    ndim = vector.ndim

    if ndim == 1:  # Handling of 1-dimensional array
        unit_vector = vector / norm
    elif ndim == 2:  # Handling of 2-dimensional array
        for i in range(0, vector.shape[1]):
            unit_vector[:, i] = vector[:, i] / norm
    else:
        log.fatal("Dimension of vector should be either 1- or 2-dimensional and not {}-dimensional.", ndim)

    return unit_vector


def _lagrange2(x, y):
    """Computes the lagrange polynomial passing through a certain set of points

    See https://en.wikipedia.org/wiki/Lagrange_polynomial and https://gist.github.com/melpomene/2482930

    Args:
        x (numpy.ndarray)      - 1-dimensional array with original x-values (e.g. time entries in [s]).
        y (numpy.ndarray)      - 1-dimensional array with original y-values (e.g. orbit position vectors).

    Returns:
        func: Lagrange interpolation function
    """

    def P(x_ip):
        total = 0
        n = len(x)
        for i in range(0, n):

            def g(i, n):
                tot_mul = 1
                for j in range(0, n):
                    if i == j:
                        continue
                    if x[i] == x[j]:
                        log.fatal(
                            "Leads to division by zero (x = {}). Identical values given in x array. For example by using Lagrange interpolation for precise orbit, check if identical observation epochs are given in SP3 file.",
                            x[i],
                        )
                    tot_mul *= (x_ip - x[j]) / float(x[i] - x[j])
                return tot_mul

            total += y[i] * g(i, n)
        return total

    return P


def _lagrange(x, y):
    """Computes the lagrange polynomial passing through a certain set of points

    Note, the performance of this routine is better than _lagrange2() or scipy.interpolate.lagrange().

    Args:
        x (numpy.ndarray)      - 1-dimensional array with original x-values (e.g. time entries in [s]).
        y (numpy.ndarray)      - 1-dimensional array with original y-values (e.g. orbit position vectors).

    Returns:
        func: Lagrange interpolation function
    """
    x_col, x_row = np.meshgrid(x, x)
    diff_x = x_col - x_row + np.eye(len(x))
    n = len(x)
    indices = np.eye(n) == 0

    def P(x_ip):
        l = np.zeros((n, x_ip.shape[0]))
        for idx, idxs in enumerate(indices):
            l[idx] = np.prod((x_ip[:, None] - x[idxs]) / diff_x[idxs, idx], axis=1)

        return np.sum(y * l.T, axis=1)

    return P


def _plot_interpolation(x, y, x_new, y_new, title=""):
    """Plot interpolated data compared to original data

    Args:
        x (numpy.ndarray)      - 1-dimensional array with original x-values (e.g. time entries in [s]).
        y (numpy.ndarray)      - n-dimensional array with original y-values (e.g. orbit position vectors).
        x_new (numpy.ndarray)  - 1-dimensional array with x-values (e.g. time entries in [s]), for which interpolation
                                 should be done.
        y_new (numpy.ndarray)  - n-dimensional array with interpolated y-values (e.g. orbit position vectors).
        title (str)            - Title of plot.
    """
    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    axes = (ax1, ax2, ax3)
    coord = ["X", "Y", "Z"]

    for idx, ax in enumerate(axes):
        ax.set_title(title + " (" + coord[idx] + " coordinate)", fontsize=12)
        ax.set_ylabel("m")
        ax.plot(x, y[:, idx], "bo", label="Original data")
        ax.plot(x_new, y_new[:, idx], "ro", label="Interpolated data")

    ax3.set_xlabel("Time")
    ax1.legend(fontsize=8, loc=1)
    f.subplots_adjust(hspace=0.3)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.show()


# def _plot_orbit_polynomial(sat, poly_coeff, x_val, y_val):
#    """Plot X, Y and Z-coordinate of given satellite orbit and the polynomial fit
#
#    Args:
#        sat         - satellite number
#        poly_coeff  - Array with coefficients of polynomial fit of SP3 orbits for each
#                      coordinate component
#        x_val       - Array with time entries
#        y_val       - 3-dim array with X, Y and Z coordinates of satellite orbit
#    """

#    f, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
#    axes = (ax1, ax2, ax3)
#    coord = ['X','Y','Z']

#    pos_fit = poly.polyval(x_val,poly_coeff).T
#    diff = pos_fit - y_val

#    for idx, ax in enumerate(axes):
#        rms = np.sqrt(np.mean(np.square(diff[:,idx])))
#        ax.set_title(coord[idx]+'-coordinate of '+sat+' orbit (RMS: '+str(rms)+' m)', fontsize=12)
#        ax.set_ylabel('m')
#        ax.plot(x_val,y_val[:,idx], 'bo',label='Original SP3 orbit data')
#        ax.plot(x_val,pos_fit[:,idx], 'ro', label='Polynominal fit of SP3 orbit data')
#        ax.plot(x_val,diff[:,idx],'go', label='Difference between fitted and original data')

#    ax3.set_xlabel('Seconds since rundate')
#    ax1.legend(fontsize=8,loc=1)
#    f.subplots_adjust(hspace=0.3)
#    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
#    plt.show()


# def triginometric_fit(x, y):
#    """Use a triginometric function for fitting given data.
#
#    Args:
#        x       - Array with x-values
#        y       - Array with y-values
#    """

#    N = 1000 # number of data points
#    t = np.linspace(0, 4*np.pi, N)
#    data = 3.0*np.sin(t+0.001) + 0.5 + np.random.randn(N) # create artificial data with noise

#    guess_mean = np.mean(data)
#    guess_std = 3*np.std(data)/(2**0.5)
#    guess_phase = 0

#    # we'll use this to plot our first estimate. This might already be good enough for you
#    data_first_guess = guess_std*np.sin(t+guess_phase) + guess_mean

#    # Define the function to optimize, in this case, we want to minimize the difference
#    # between the actual data and our "guessed" parameters
#    optimize_func = lambda x: x[0]*np.sin(t+x[1]) + x[2] - data
#    est_std, est_phase, est_mean = leastsq(optimize_func, [guess_std, guess_phase, guess_mean])[0]

#    # recreate the fitted curve using the optimized parameters
#    data_fit = est_std*np.sin(t+est_phase) + est_mean

#    plt.plot(data, '.')
#    plt.plot(data_fit, label='after fitting')
#    plt.plot(data_first_guess, label='first guess')
#    plt.legend()
#    plt.show()


def rotate_x(angle):
    """Rotation matrix around X-axis

    Positive (counterclockwise) rotation of the X-axis as viewed from the positive end of the rotation axis towards
    the origin.

    Args:
        angle (float64):    Rotation angle in [rad]

    Return:
        numpy.ndarray:      Rotation matrix
    """
    log.dev("lib.mathp.rotate_x is deprecated. Use lib.rotation.R1 instead.")

    cosA = np.cos(angle)
    sinA = np.sin(angle)
    R = np.array([[1, 0, 0], [0, cosA, sinA], [0, -sinA, cosA]])
    return R


def rotate_y(angle):
    """Rotation matrix around Y-axis

    Positive (counterclockwise) rotation of the Y-axis as viewed from the positive end of the rotation axis towards
    the origin.

    Args:
        angle (float64):    Rotation angle in [rad]

    Return:
        numpy.ndarray:      Rotation matrix
    """
    log.dev("lib.mathp.rotate_y is deprecated. Use lib.rotation.R2 instead.")

    cosA = np.cos(angle)
    sinA = np.sin(angle)
    R = np.array([[cosA, 0, -sinA], [0, 1, 0], [sinA, 0, cosA]])
    return R


def rotate_z(angle):
    """Rotation matrix around Z-axis

    Positive (counterclockwise) rotation of the Z-axis as viewed from the positive end of the rotation axis towards
    the origin.

    Args:
        angle (float64):    Rotation angle in [rad]

    Return:
        numpy.ndarray:      Rotation matrix
    """
    log.dev("lib.mathp.rotate_z is deprecated. Use lib.rotation.R3 instead.")

    cosA = np.cos(angle)
    sinA = np.sin(angle)
    R = np.array([[cosA, sinA, 0], [-sinA, cosA, 0], [0, 0, 1]])
    return R
