# cython: profile=True
# cython: linetrace=True
# distutils: define_macros=CYTHON_TRACE_NOGIL=1
"""Framework for calculating satellite orbit models

Description:
------------

Each orbit model should be defined in a separate .py-file. The function inside the .py-file that should be called need
to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    def gravity_earth(...):
        ...

The decorated function will be called with several parameters::

    TODO: Document parameters


References:
-----------

[1] Montenbruck, Oliver and Gill, Eberhard: Satellite Orbits, Springer Verlag, 2000.
"""

# Standard library imports
import os
import sys
cimport libc.math
import math
from collections import OrderedDict

# External library imports
import numpy as np

# Midgard imports
from midgard.dev.timer import Timer
from midgard.dev import plugins
from midgard.dev import log
from midgard.math import interpolation

# Where imports
from where import apriori
from where.lib import config
from where.lib.time import Time
from where.lib.time import TimeDelta
from where import integrators

_INITIAL_POSVEL = dict()
_TRANSITION_MATRIX = dict()
_PARAMETERS = dict()
_MEASUREMENT_PARAMETERS = dict()
_SENSITIVITY_MATRIX = dict()
_MEASUREMENT_MATRIX = dict()
_MODELS = dict()


def calculate(rundate, sat_name, obs_sec, range_bias_stations, time_bias_stations, return_full_table=False):
    """Calculate the orbit for a satellite at the given observation epochs

    Solve the differential equation of motion of the satellite and the differential equation of the state transition
    matrix of the satellite simultaneously.  With 6 variables, 3 for position and 3 for velocity, in the state vector
    the state transition matrix becomes 6x6, so in total we need to solve for 42 variables.

    The actual equations are defined in the gravity_field-specific code inside the gravity-package. Although general
    corrections, e.g. due to the gravitational forces of the Moon and the Sun are defined in this module.

    The return values from the function depends on the return_full_table-parameter. If return_full_table is False
    (default), then only arrays of satellite position and velocity are returned with values at the provided observation
    epochs. If return_full_table is True then arrays of satellite position and velocity is returned in addition to an
    array containing a time grid in seconds since rundate. Values are provided for the full time grid which will depend
    on the orbit_step_length and the observation epochs.

    Args:
        rundate:           The model run date.
        sat_name:          Name of the satellite.
        obs_sec:           List of observation epochs in seconds after rundate.
        return_full_table: Arrays are for observation epochs or full time grid.

    Returns:
        Arrays of calculated positions and velocities.
    """
    # Read configuration settings
    integrate_method = config.tech.integrate_method.str

    # Set initial state of satellite orbit: position and velocity, and
    # of the state transition matrix, the 6x6 identity matrix
    if sat_name not in _INITIAL_POSVEL:
        satellite = apriori.get_satellite(*sat_name.split("-"))
        _INITIAL_POSVEL[sat_name] = satellite.initial_posvel(rundate)

    # Set parameters to be estimated:
    force_parameter_names = config.tech.force_parameters.list
    measurement_parameter_names = config.tech.measurement_parameters.list
    set_parameters(sat_name, force_parameter_names, rundate)
    set_measurement_parameters(sat_name, measurement_parameter_names, rundate, range_bias_stations, time_bias_stations)
    num_param = len(_PARAMETERS[sat_name])
    end_time = max(obs_sec)
    log.info(
        f"Calculating orbit of {sat_name} from {rundate.strftime(config.FMT_datetime)} to "
        f"{(Time(rundate, scale='utc') + TimeDelta(end_time, format='sec'))}"
    )

    step_length = config.tech.orbit_step_length.float
    arc_length = config.tech.arc_length.int
    grid_step = step_length if step_length else 10.0

    initial_state = np.hstack((_INITIAL_POSVEL[sat_name], np.identity(6).reshape(-1), np.zeros(6 * num_param)))

    # Set up time grid
    c = 14 # Corresponds to constant in integrator cowell
    time_grid2 = np.arange(-c * grid_step, math.ceil(end_time / grid_step) * grid_step + grid_step, grid_step)

    # Solve the equation of motion for the state vector and the state transition matrix simultaneously. In total, 42
    # first order differential equations are solved by the integrator.
    integrand = construct_forces(rundate, sat_name, arc_length, num_param, time_grid2, c)

    log.info(f"Integrating orbit to {grid_step}-second time grid")
    log.info("Initial position: {:0.8f} {:0.8f} {:0.8f}".format(*initial_state[0: 3]))
    log.info("Initial velocity: {:0.8f} {:0.8f} {:0.8f}".format(*initial_state[3: 6]))

    with Timer('Finish integrating orbit in'):
        state, time_grid = integrators.call(integrate_method, integrand=integrand, initial_state=initial_state,
                                            grid_step=grid_step, end_time=end_time)
    if not set(time_grid).issubset(time_grid2):
        log.fatal("Something wrong with the time grid")

    sat_pos = state[:, :3]
    sat_vel = state[:, 3: 6]
    trans_mat = interpolation.interpolate(np.array(time_grid), state[:, 6: 42], obs_sec,
                                                     kind="interpolated_univariate_spline")
    # Store the state transition matrix at observation epochs for
    # updating the initial state in later iterations.
    _TRANSITION_MATRIX[sat_name] = trans_mat.reshape(-1, 6, 6)
    # Store the sensitivity matrix
    if num_param > 0:
        sens_mat = interpolation.interpolate(np.array(time_grid), state[:, 42:], obs_sec,
                                                        kind="interpolated_univariate_spline")
        _SENSITIVITY_MATRIX[sat_name] = sens_mat.reshape(-1, 6, num_param)

    if not return_full_table:
        sat_pos = interpolation.interpolate(np.array(time_grid), sat_pos, obs_sec,
                                                                      kind="interpolated_univariate_spline")
        sat_vel = interpolation.interpolate(np.array(time_grid), sat_vel, obs_sec,
                                                                      kind="interpolated_univariate_spline")
        time_grid = obs_sec

    return sat_pos, sat_vel, time_grid
    #return sat_pos[:-c, :], sat_vel[:-c, :], time_grid[c:]


cdef construct_forces(rundate, sat_name, int arc_length, int num_param, double[:] time_grid, int c):
    """Construct the forces needed to find the satellite orbit

    The forces_func below will be the \f$ \vec F \f$ in the Newtonian differential equation

    \f[ \ddot{\vec r} = \vec F(t, \vec r, \dot{\vec r}) / m , \f]

    where \f$ \vec r \f$ is the position of the satellite. That is, \f$ F \f$ describes the forces acting on the
    satellite. To solve the equation, we rewrite it to a first order differential equation by letting \f$ \vec y =
    (\vec r, \dot{\vec r}) \f$ be the vector combining both position and velocity. Then we can write

    \f[ \dot{\vec y} = \vec F(t, \vec y) / m . \f]

    The definition of \f$ \vec F \f$ has been expanded to include the identity \f$ \dot{\vec r} = \dot{\vec r} \f$, so
    that we have a system of six differential equations (with the notation \f$ \vec r = (r_x, r_y, r_z) \f$ and \f$
    \dot{\vec r} = (v_x, v_y, v_z) \f$),

    \f[ \left( \begin{array}{c} \dot{r_x} \\ \dot{r_y} \\ \dot{r_z} \\
      \dot{v_x} \\ \dot{v_y} \\ \dot{v_z} \end{array} \right) = \left(
      \begin{array}{c} v_x \\ v_y \\ v_z \\ F_x(t, \vec y) / m \\
      F_y(t, \vec y) / m \\ F_z(t, \vec y) / m \end{array} \right) . \f]

    Args:
        rundate:                         The rundate of the model.
        sat_name:                        Name of the satellite.
        arc_length:                      Number of days to integrate over.
        num_param:                       Number of parameters to be estimated
    Returns:
        A function representing the \f$ F \f$-integrand.

    """
    cdef int i, j, idx
    cdef list bodies
    cdef int num_bodies
    vars_ = dict(rundate=rundate, sat_name=sat_name, force_parameters=_PARAMETERS[sat_name], num_param=num_param)

    epochs = Time([Time(rundate, scale="utc") + TimeDelta(t, format="sec") for t in time_grid])
    cdef double[:, :, :] gcrs2itrs = epochs.gcrs2itrs

    # Ephemerides of sun, moon and planets
    bodies = config.tech.gravity_bodies.bodies.list
    num_bodies = len(bodies)
    cdef double[:, :, :] body_pos_gcrs = np.zeros((num_bodies, len(epochs), 3))
    cdef double[:, :, :] body_pos_itrs = np.zeros((num_bodies, len(epochs), 3))

    eph = apriori.get("ephemerides")
    for idx in range(0, num_bodies):
        pos_gcrs = eph.pos_gcrs(bodies[idx], time=epochs)
        pos_itrs = eph.pos_itrs(bodies[idx], time=epochs)
        for i in range(0, len(epochs)):
            for j in range(0, 3):
                body_pos_gcrs[idx, i, j] = pos_gcrs[i, j]
                body_pos_itrs[idx, i, j] = pos_itrs[i, j]

    read_models(rundate, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs)

    def forces_func(double[:] state, int n):
        """The function that is integrated to find the satellite orbit

        Args:
            state:    42 + 6 * num_param floats, position, velocity, transition matrix and sensitivity matrix
            n:        Integer, where the current step of the integrator is n + c. Hence n = 0 corresponds to rundate.
        Returns:
            Vector of 42  + 6 * num_param floats, right hand side of differential equation.
        """
        # Update variables
        vars_.update(dict(sat_pos_gcrs=state[: 3], sat_vel_gcrs=state[3: 6], sat_pos_itrs=np.dot(gcrs2itrs[n + c],
                          state[: 3]), sat_vel_itrs=np.dot(gcrs2itrs[n + c], state[3: 6]), current_step=n + c))

        # Calculate the right hand side of the variational equation for the acceleration, transition matrix and
        # sensitivity matrix for each force
        cdef double[:, :] trans_matrix = np.zeros((6, 6))
        cdef int i, j, k, l
        cdef double[:,:] tr, s
        cdef double[:] sat_acc_gcrs

        trans_matrix[0, 3] = 1
        trans_matrix[1, 4] = 1
        trans_matrix[2, 5] = 1
        # trans_matrix[:3, 3:] = np.eye(3)

        cdef double[:, :] sens_matrix = np.zeros((6, num_param))

        sat_acc_gcrs, tr, s = calculate_epoch(**vars_)
        for i in range(0, 3):
            for j in range(6):
                trans_matrix[3 + i, j] = tr[i, j]
        for i in range(0, 3):
            for j in range(num_param):
                sens_matrix[3 + i, j] = s[i, j]
        # phi_matrix = state[6:42].reshape(6, 6)
        # s_matrix = state[42:].reshape(6, num_param)

        cdef double[:, :] phi_matrix = np.zeros((6, 6))
        cdef double[:, :] s_matrix = np.zeros((6, num_param))
        k = 0
        for i in range(0, 6):
            for j in range(0, 6):
                phi_matrix[i, j] = state[6 + k]
                k = k + 1
        k = 0
        for i in range(0, 6):
            for j in range(0, num_param):
                s_matrix[i, j] = state[42 + k]
                k = k + 1

        cdef double[:] vec = np.zeros(42 + 6 * num_param)
        for i in range(3):
            vec[i] = vars_["sat_vel_gcrs"][i]
            vec[3 + i] = sat_acc_gcrs[i]
        l = 6
        # np.dot(trans_matrix, phi_matrix).reshape(-1)
        for i in range(6):
            for j in range(6):
                for k in range(6):
                    vec[l] += trans_matrix[i, k] * phi_matrix[k, j]
                l = l + 1

        # np.dot(trans_matrix, s_matrix).reshape(-1) + sens_matrix.reshape(-1)
        for i in range(6):
            for j in range(num_param):
                vec[l] = sens_matrix[i, j]
                for k in range(6):
                    vec[l] += trans_matrix[i, k] * s_matrix[k, j]
                l = l + 1
        # Return before cythonizing:
        # return velocity, sat_acc_gcrs, trans_matrix * phi_matrix, trans_matrix * s_matrix + sens_matrix
        return vec

    return forces_func


cdef set_parameters(sat_name, parameter_names, rundate):

    orbit_models = config.tech.orbit_models.list

    if "empirical" in orbit_models and "empirical" not in parameter_names:
        log.fatal("Empirical is included in orbit_models, not in force_parameters")
    if "empirical" in parameter_names:
        empirical_parameters = config.tech.empirical.empirical_parameters.list
        parameter_names.remove("empirical")
        parameter_names = parameter_names + empirical_parameters

    if sat_name not in _PARAMETERS:
        satellite = apriori.get_satellite(*sat_name.split("-"))
        _PARAMETERS[sat_name] = OrderedDict()

        for i in range(0, len(parameter_names)):
            # Test if parameters given on config file are valid
            # and set default values for parameters
            if parameter_names[i] == "drag_coefficient":
                if "drag" in orbit_models:
                    _PARAMETERS[sat_name]["drag_coefficient"] = satellite.drag_product
                    continue
                else:
                    log.warn("Could not estimate drag coefficient")
                    log.fatal("Drag is not included in orbit models")
            elif parameter_names[i] == "radiation_pressure_coefficient":
                if "solar_radiation_pressure" in orbit_models:
                    _PARAMETERS[sat_name]["radiation_pressure_coefficient"] = satellite.radiation_pressure_coefficient
                    continue
                else:
                    log.warn("Could not estimate radiation pressure coefficient")
                    log.fatal("Solar radiation pressure is not included in orbit models")
            elif parameter_names[i].startswith("a"):
                _PARAMETERS[sat_name][parameter_names[i]] = 0
            elif parameter_names[i] == "c20":
                if "gravity_earth" in orbit_models:
                    gravity_field = config.tech.gravity_field.str
                    gravity_coeffs = apriori.get(
                        "gravity", gravity_field=gravity_field, truncation_level=2, rundate=rundate
                    )
                    _PARAMETERS[sat_name]["c20"] = gravity_coeffs["C"][2, 0]
                    log.info("in order to get initial value for gravity coefficient")
                    log.info(f"C20 = {_PARAMETERS[sat_name]['c20']}")
                    continue
                else:
                    log.warn("Could not estimate c20")
                    log.fatal("gravity_earth is not included in orbit models")
            else:
                log.fatal(f"Not a valid parameter name {parameter_names[i]!r}")


cdef set_measurement_parameters(sat_name, parameter_names, rundate, range_bias_stations, time_bias_stations):
    # TODO: Add station in dict
    if sat_name not in _MEASUREMENT_PARAMETERS:
        satellite = apriori.get_satellite(*sat_name.split("-"))
        _MEASUREMENT_PARAMETERS[sat_name] = OrderedDict()
        for i in range(0, len(parameter_names)):
            # Test if parameters given on config file are valid
            # and set default values for parameters
            if parameter_names[i] == "range_bias":
                for j in range(0, len(range_bias_stations)):
                    _MEASUREMENT_PARAMETERS[sat_name][f"range_bias_{range_bias_stations[j]}"] = 0
            elif parameter_names[i] == "time_bias":
                for j in range(0, len(time_bias_stations)):
                    _MEASUREMENT_PARAMETERS[sat_name][f"time_bias_{time_bias_stations[j]}"] = 0
            else:
                log.fatal(f"Not a valid parameter name {parameter_names[i]!r}")



def update_orbit(sat_name, sta_pos_gcrs, sat_pos_gcrs, sat_vel, residual, weight, stations):
    """
    Compute correction to the initial value of position
    and velocity of the satellite in meters and meters/sec

    Args:
        sat_name:     Name of satellite
        sta_pos_gcrs: Station position in the gcrs system
        sat_pos_gcrs: Satellite position in the gcrs system
        residual:     Observation residual
        weight:       Weight of the observation
        stations:     List of stations, one for each observations

    References:
         Montenbruck and Gill [1], Section 8.1.1.
    """
    num_obs = len(weight)
    sta_sat_vector = sat_pos_gcrs - sta_pos_gcrs
    unit_vector = (sta_sat_vector / np.linalg.norm(sta_sat_vector, axis=1)[:, None])
    unit_vector_ext = np.hstack((unit_vector, np.zeros(unit_vector.shape)))
    weight_matrix = np.diag(weight)
    
    if sat_name in _SENSITIVITY_MATRIX:
        H = np.hstack((np.sum(unit_vector_ext[:, :, None] * _TRANSITION_MATRIX[sat_name], axis=1),
                       np.sum(unit_vector_ext[:, :, None] * _SENSITIVITY_MATRIX[sat_name], axis=1)))
    else:
        H = np.sum(unit_vector_ext[:, :, None] * _TRANSITION_MATRIX[sat_name], axis=1)

    _MEASUREMENT_MATRIX[sat_name] = np.zeros((len(stations), len(_MEASUREMENT_PARAMETERS[sat_name])))
    for i in range(0, len(stations)):
        for j, (k, val) in enumerate(_MEASUREMENT_PARAMETERS[sat_name].items()):
            if f"range_bias_{stations[i]}" == k:
                _MEASUREMENT_MATRIX[sat_name][i, j] = 1
            elif f"time_bias_{stations[i]}" == k:
                # TODO! Not sure if this is the correct way of estimating partials of 
                # measurement with respect to time bias
                _MEASUREMENT_MATRIX[sat_name][i, j] = np.sum(unit_vector * sat_vel, axis=1)[i]
    H = np.hstack((H, _MEASUREMENT_MATRIX[sat_name]))
    
    # Equation 8.15 in Montenbruck and Gill:
    N = np.linalg.pinv(np.dot(H.T, weight_matrix @ H))
    NH = np.dot(N, H.T)
    _INITIAL_POSVEL[sat_name] += np.dot(NH, weight_matrix @ residual)[: 6]

    log.info(f"Estimated initial position {_INITIAL_POSVEL[sat_name][0:3]}")
    log.info(f"Estimated initial velocity {_INITIAL_POSVEL[sat_name][3:6]}")

    if sat_name in _SENSITIVITY_MATRIX:
        keys = list(_PARAMETERS[sat_name].keys())
        keys_measurement = list(_MEASUREMENT_PARAMETERS[sat_name].keys())
        for i in range(0, len(keys)):
            _PARAMETERS[sat_name][keys[i]] += np.dot(NH, weight_matrix @ residual)[6:][i]
            log.info(f"Estimate of {keys[i]} is {_PARAMETERS[sat_name][keys[i]]}:")
        for i in range(0, len(keys_measurement)): 
            _MEASUREMENT_PARAMETERS[sat_name][keys_measurement[i]] += np.dot(NH, weight_matrix @ residual)[6 + len(keys):][i]
            log.info(f"Estimate of {keys_measurement[i]} is {_MEASUREMENT_PARAMETERS[sat_name][keys_measurement[i]]}:")

    estimated_range_bias = np.zeros(num_obs)
    estimated_time_bias = np.zeros(num_obs)
    for i in range(0, len(stations)):
        for (k, val) in _MEASUREMENT_PARAMETERS[sat_name].items():
            if f"range_bias_{stations[i]}" == k:
                estimated_range_bias[i] = val
            if f"time_bias_{stations[i]}" == k:
                estimated_time_bias[i] = val
    return estimated_range_bias, estimated_time_bias

def read_models(rundate, sat_name, time_grid, epochs, body_pos_gcrs, body_pos_itrs, bodies, gcrs2itrs):
    """Read Cython or Python orbit models
    Args:
        rundate:           Time of integration start
        sat_name:          Name of satellite
        time_grid:         Table of times in seconds since rundate, in utc.
        epochs:            time_grid converted to Time objects, in utc.
        body_pos_gcrs:     The positions of the bodies in the solar system in GCRS.
        body_pos_itrs:     The positions of the bodies in the solar system in ITRS.
        bodies:            List of bodies
        gcrs2itrs:         List of transformation matrices, one for each time in epochs.
    """
    orbit_models = config.tech.orbit_models.list

    log.info("Using Cython implementation of orbit models")

    import importlib
    package = __name__.rsplit(".", maxsplit=1)[0]
    modules = [m[:-4] for m in os.listdir(os.path.dirname(__file__))
               if (m.endswith(".pyx") and not m.startswith("_"))]
    for module_name in modules:
        try:
            module = importlib.import_module("{}.{}".format(package, module_name))
            entry_points = module.register_entry_point()
            if "setup" in entry_points and module_name in orbit_models:
                entry_points["setup"](
                    rundate,
                    _PARAMETERS[sat_name],
                    sat_name,
                    time_grid,
                    epochs,
                    body_pos_gcrs,
                    body_pos_itrs,
                    bodies,
                    gcrs2itrs,
                )
            _MODELS[module_name] = entry_points["call"]
        except ImportError:
            log.warn(f"Did not find compiled Cython module '{package}.{module_name}'")

    # Compare with orbit_models
    for missing_model in (set(orbit_models) - set(_MODELS.keys())):
        log.warn(f"Unknown model {missing_model!r}")

    for not_included_model in (set(_MODELS.keys()) - set(orbit_models)):
        del _MODELS[not_included_model]

    log.info(f"Integrating with orbit models {', '.join(_MODELS)}")


def calculate_epoch(**kwargs):
    """Call all orbit models

    The list of models to apply is taken from the config file of the
    given technique.

    Args:
        kwargs:  Variables that are passed on to the individual models.

    Returns:
        Acceleration (3-vector) and lower 3x6 transition matrix.
    """
    sat_name = kwargs["sat_name"]   # TODO: sat_name is not used
    num_param = kwargs["num_param"]
    # Call each orbit model
    acceleration = np.zeros(3)
    transition_matrix = np.zeros((3, 6))
    sensitivity_matrix = np.zeros((3, num_param))
    for model, model_func in _MODELS.items():
        acc_part, trans_part, sens_part = model_func(**kwargs)
        acceleration += acc_part
        transition_matrix += trans_part
        sensitivity_matrix += sens_part

    return acceleration, transition_matrix, sensitivity_matrix


    

