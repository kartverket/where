"""Correct clocks toward a reference clock

Description:
------------

Corrects the station clocks by using a quadratic (3 terms) clock polynomial applied to each station except one which is
taken to be the reference station. Uses least squares to minimize the residuals.

If a station has one or multiple clock breaks, separate polynomials are estimated for the data before and after the
clock break.




"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.data.time import Time, TimeDelta

# Name of model
MODEL = __name__.split(".")[-1]


@plugins.register
def clock_correction(dset):
    """Estimate clock polynomial
    """
    # Take previous clock corrections into account
    try:
        output = dset.vlbi_clock
    except AttributeError:
        output = np.zeros(dset.num_obs)

    # Read order of clock polynomial from config file
    terms = 1 + config.tech.get("order_of_polynomial", section=MODEL, default=2).int

    # Read clock breaks from session config, only split on commas (and commas followed by whitespace)
    clock_breaks = config.tech.get("clock_breaks", section=MODEL, default="").as_list(split_re=", *")
    stations, time_intervals = parse_clock_breaks(dset, clock_breaks)

    # Read reference clock from edit file and store in dataset
    ref_clock_str = config.tech.get("reference_clock", section=MODEL, default="").str
    ref_clock = parse_reference_clock(stations, ref_clock_str)
    dset.meta["ref_clock"] = ref_clock

    # Remove reference clock from list of clocks to be estimated
    idx = stations.index(ref_clock)
    del stations[idx]
    del time_intervals[idx]

    # Number of clock polynomial coefficients
    num_coefficients = len(stations) * terms
    param_names = [
        sta + " clk_a" + str(t) + " " + time_intervals[i][0].utc.iso + " " + time_intervals[i][1].utc.iso
        for i, sta in enumerate(stations)
        for t in range(terms)
    ]
    dset.meta["num_clock_coeff"] = num_coefficients

    # Set up matrices for estimation
    A = np.zeros((dset.num_obs, num_coefficients, 1))

    # Time coefficients, used when setting up A
    t = dset.time.utc.mjd - dset.time.utc[0].mjd
    poly = np.array([t ** n for n in range(terms)]).T

    # Set up the A matrix with time coefficients
    for idx, (station, (t_start, t_end)) in enumerate(zip(stations, time_intervals)):
        filter_time = np.logical_and(t_start.utc.mjd <= dset.time.utc.mjd, dset.time.utc.mjd < t_end.utc.mjd)
        filter_1 = np.logical_and(dset.filter(station_1=station), filter_time)
        A[filter_1, idx * terms : (idx + 1) * terms, 0] = poly[filter_1]
        filter_2 = np.logical_and(dset.filter(station_2=station), filter_time)
        A[filter_2, idx * terms : (idx + 1) * terms, 0] = -poly[filter_2]

    # Calculate normal matrix N and the moment vector U
    U = np.sum(A @ dset.residual[:, None, None], axis=0)
    N = np.sum(A @ A.transpose(0, 2, 1), axis=0)

    # Invert the normal matrix to find corrections, only the non-zero part of the matrix is inverted
    idx = np.logical_not(U == 0)[:, 0]
    X = np.zeros((num_coefficients, 1))
    det = np.linalg.det(N[idx, :][:, idx])
    threshold = 1e-12
    if np.abs(det) < threshold:
        # TODO: what is a good threshold value?
        rank = np.linalg.matrix_rank(N[idx, :][:, idx])
        log.warn(f"Determinant of normal matrix in clock correction is close to zero ({det})")
        log.info(f"Normal matrix shape = {N.shape}, normal matrix rank = {rank}")
        _, R = np.linalg.qr(N[idx, :][:, idx])
        for i, row in enumerate(R):
            if np.max(np.abs(row)) < threshold * 10 ** 3:
                log.error(f"{param_names[i]} linearly dependent (max_row = {np.max(np.abs(row))})")
    try:
        X[idx] = np.linalg.inv(N[idx, :][:, idx]) @ U[idx]
    except np.linalg.LinAlgError:
        log.fatal(f"Singular matrix in {MODEL}")

    # Calculate final corrections
    output += (A.transpose(0, 2, 1) @ X)[:, 0, 0]
    return output


def parse_clock_breaks(dset, clock_breaks):
    """Parses the clock breaks string from the edit file

    Examples:
        > parse_clock_breaks(dset, '')
        (OrderedDict(), 0)
        > parse_clock_breaks(dset, 'SVETLOE 2015/01/23 05:36:00,
                                    SVETLOE 2015/01/23 16:53:00,
                                    SVETLOE 2015/01/23 12:30:00')
        (OrderedDict([(5, [106560.0, 131400.0, 147180.0])]), 3)

    Args:
        dset:                A Dataset containing model data.
        clock_breaks_str:    A string with clock break information

    Returns:
        OrderedDict with clock breaks and total number of clock breaks
     """
    # Parse clock breaks from file and store in the station_breaks dictionary
    station_breaks = {
        s: [min(dset.time.utc), max(dset.time.utc) + TimeDelta(1, fmt="seconds", scale="utc")]
        for s in dset.unique("station")
    }
    if clock_breaks:
        log.info(f"Applying clock breaks: {', '.join(clock_breaks)}")

    for cb in clock_breaks:
        # Station names may contain spaces
        cb = cb.split()
        cb_date = cb[-2:]
        cb_station = " ".join(cb[:-2])
        # cb_station, *cb_date = cb.split()
        cb_time = Time(" ".join(cb_date), scale="utc", fmt="iso")
        if cb_station not in station_breaks:
            log.warn(
                f"Station {cb_station} with clock break unknown. Available options are {', '.join(station_breaks)}"
            )
            continue
        station_breaks[cb_station].append(cb_time)
        dset.meta.add_event(cb_time, "clock_break", cb_station)

    # Convert the station_breaks dict to lists of (station, (time_start, time_end))-tuples
    stations = list()
    time_intervals = list()
    for station in sorted(station_breaks.keys(), key=lambda s: (len(station_breaks[s]), s), reverse=True):
        station_times = sorted(station_breaks[station])
        for t_start, t_end in zip(station_times[:-1], station_times[1:]):
            stations.append(station)
            time_intervals.append((t_start, t_end))

    return stations, time_intervals


def parse_reference_clock(stations, ref_clock_str):
    """Parses the reference clock string from the edit file

    Args:
        dset:              A Dataset containing model data
        ref_clock_str:     IVS name of reference clock station

    Returns:
        String: IVS name of reference clock station in Dataset
    """
    if ref_clock_str not in stations:
        if ref_clock_str:
            log.warn(f"Reference clock {ref_clock_str!r} unknown. Available options are {', '.join(stations)}")

        # Pick last station as default
        ref_clock_str = stations[-1]
    log.info(f"Reference clock is {ref_clock_str!r}")
    return ref_clock_str
