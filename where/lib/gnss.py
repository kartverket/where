#!/usr/bin/env python3
"""Where library module including functions for GNSS modeling

Example:
--------

    from where.lib import gnss
    ...

Description:
------------

This module will provide functions for GNSS modeling.


TODO: How to move routines to Midgard?
========================================
check_satellite_eclipse(dset)
findsun(time)                                  -> Midgard: planetary_motion?
gsdtime_sun(time)                              -> Midgard: planetary_motion?
get_earth_rotation(dset)                       -> Midgard: PosVel (see function)
get_code_observation(dset)                     -> Midgard: gnss
get_flight_time(dset)                          -> Midgard: PosVel (see function)
obstype_to_freq(sys, obstype)                  -> Midgard: gnss
get_initial_flight_time(dset, sat_clock_corr=None, rel_clock_corr=None)  -> Midgard: gnss
get_line_of_sight(dset)                        -> Midgard: Position library
get_rinex_file_version(file_key, file_vars)    -> Is that needed in the future?
gpssec2jd(wwww, sec)                           -> in time.gps_ws
Example:
    from where.data import time
    t = time.Time([2000, 2000, 2000, 2004], [0, 100000, 200000, 86400], fmt="gps_ws", scale="gps")
    t.gps_ws
jd2gps(jd)                                     -> in time
linear_combination(dset)                       -> Midgard: gnss
llh2xyz(lat, lon, h)                           -> midgard.math.transformation.llh2trs
plot_skyplot(dset)                             -> Midgard: plot

Should we more specific in using arguments, that instead of using 'dset'? -> Maybe

"""
# Standard library imports
from typing import Tuple

# External library imports
import numpy as np
import matplotlib.pyplot as plt

# Migard imports
from midgard.collections import enums
from midgard.math.constant import constant

# Where imports
from where.lib import config
from where.lib import log
from where.lib import mathp
from where.lib import rotation
from where.data.time import TimeDelta


def check_satellite_eclipse(dset):
    """Check if a satellite is an eclipse

    TODO: Check if a better algorithm exists (e.g. based on beta angle).

    Args:
        dset(Dataset):    Model data
    """
    cos_gamma = np.einsum(
        "ij,ij->i", mathp.unit_vector(dset.sat_posvel.itrs_pos), dset.sat_posvel.itrs_pos_sun
    )  # TODO:  dot product -> better solution dot() function in mathp
    h = np.linalg.norm(dset.sat_posvel.itrs_pos, axis=1) * np.sqrt(1.0 - cos_gamma ** 2)

    satellites_in_eclipse = list()
    for satellite in dset.unique("satellite"):
        idx = dset.filter(satellite=satellite)
        satellite_eclipse = np.logical_and(cos_gamma[idx] < 0, h[idx] < constant.a)
        if np.any(satellite_eclipse == True):
            satellites_in_eclipse.append(satellite)

    return satellites_in_eclipse


def findsun(time):
    """Obtains the position vector of the Sun in relation to Earth (in ECEF).

    This routine is a reimplementation of routine findSun() in model.c of gLAB 3.0.0 software.

    Args:
        time(Time):    Time object

    Returns:
        numpy.ndarray:  Sun position vector given in ECEF [m]
    """
    AU = 1.495_978_70e8
    gstr, slong, sra, sdec = gsdtime_sun(time)

    sun_pos_x = np.cos(np.deg2rad(sdec)) * np.cos(np.deg2rad(sra)) * AU
    sun_pos_y = np.cos(np.deg2rad(sdec)) * np.sin(np.deg2rad(sra)) * AU
    sun_pos_z = np.sin(np.deg2rad(sdec)) * AU
    sun_pos_eci = np.vstack((sun_pos_x, sun_pos_y, sun_pos_z)).T

    # Rotate from inertial to non inertial system (ECI to ECEF)
    sun_pos_ecef = (rotation.R3(np.deg2rad(gstr)) @ sun_pos_eci.T)[:, :, 0]  # remove 1 dimension

    return sun_pos_ecef


def gsdtime_sun(time):
    """Get position of the sun (low-precision)

    This routine is a reimplementation of routine GSDtime_sun() in model.c of gLAB 3.0.0 software.

    Args:
        time(Time):    Time object

    Returns:
        tuple:  with following entries

    =============== =============== ==================================================================================
     Elements        Type            Description
    =============== =============== ==================================================================================
     gstr            numpy.ndarray   GMST0 (to go from ECEF to inertial) [deg]
     slong           numpy.ndarray   Sun longitude [deg]
     sra             numpy.ndarray   Sun right Ascension [deg]
     sdec            numpy.ndarray   Sun declination in [deg]
    =============== =============== ==================================================================================
    """
    jd = time.mjd_int - 15019.5
    frac = time.jd_frac
    vl = np.mod(279.696_678 + 0.985_647_335_4 * jd, 360)
    gstr = np.mod(279.690_983 + 0.985_647_335_4 * jd + 360 * frac + 180, 360)
    g = np.deg2rad(np.mod(358.475_845 + 0.985_600_267 * jd, 360))

    slong = vl + (1.91946 - 0.004_789 * jd / 36525) * np.sin(g) + 0.020_094 * np.sin(2 * g)
    obliq = np.deg2rad(23.45229 - 0.013_012_5 * jd / 36525)

    slp = np.deg2rad(slong - 0.005_686)
    sind = np.sin(obliq) * np.sin(slp)
    cosd = np.sqrt(1 - sind * sind)
    sdec = np.rad2deg(np.arctan2(sind, cosd))

    sra = 180 - np.rad2deg(np.arctan2(sind / cosd / np.tan(obliq), -np.cos(slp) / cosd))

    return gstr, slong, sra, sdec


# TODO: pv.trs.observed - pv.trs # calculate property 'observed' = rotation.R3(rotation_angle[idx]).dot(dset.sat_posvel.itrs_pos[idx])
# def get_earth_rotation(posvel: PositionVelocityArray, flight_time: np.ndarray):
def get_earth_rotation(dset):
    """Get corrections for satellite position and velocity by Earth rotation

    In a Earth-fixed reference system the Earth's rotation has to be applied, which accounts for time effect of Earth
    rotation during the signal propagates from the satellite to the receiver. Eq. 5.11 in :cite:`subirana2013` is used
    for correcting the satellite position and velocity in the Dataset field 'sat_posvel' about the Earth's rotation
    effect.

    Args:
        dset(Dataset):    Model data

    Returns:
        tuple:    with following entries

    =============== =============== ==================================================================================
     Elements        Type            Description
    =============== =============== ==================================================================================
     sat_pos         numpy.ndarray   Satellite position vector corrections in ITRS and [m]
     vel_pos         numpy.ndarray   Satellite velocity corrections in ITRS and [m/s]
    =============== =============== ==================================================================================
    """
    sat_pos = np.zeros((dset.num_obs, 3))
    sat_vel = np.zeros((dset.num_obs, 3))

    flight_time = get_flight_time(dset)
    rotation_angle = flight_time * constant.omega

    sat_pos = (
        rotation.R3(rotation_angle) @ dset.sat_posvel.trs.pos.val[:, :, None] - dset.sat_posvel.trs.pos.val[:, :, None]
    )
    sat_vel = (
        rotation.R3(rotation_angle) @ dset.sat_posvel.trs.vel.val[:, :, None] - dset.sat_posvel.trs.vel.val[:, :, None]
    )

    return sat_pos[:, :, 0], sat_vel[:, :, 0]


def get_code_observation(dset):
    """Get pseudo-range (code) observations depending on given observation types

    The first element of the observation type variable `dset.meta['obstypes'][sys]` is selected as observation for
    single frequency solution. The order of the observation type variable `dset.meta['obstypes'][sys]` depends on
    the priority list given in the configuration file and the given observations.

    The ionospheric-free linear combination is applied for dual frequency solution.

    Args:
        dset:    Dataset

    Returns:
        numpy.ndarray:  Pseudo-range (code) observation choosen depending on priority list and for dual frequency
                        solution given as ionospheric-free linear combination
    """
    freq_type = config.tech.freq_type.str
    code_obs = np.zeros(dset.num_obs)

    if freq_type == "single":

        for sys in dset.unique("system"):
            idx = dset.filter(system=sys)
            obstype = dset.meta["obstypes"][sys][0]
            code_obs = dset.obs[obstype][idx]

    elif freq_type == "dual":
        code_obs, _ = linear_combination("ionosphere-free", dset)
    else:
        log.fatal(
            "Configuration option 'freq_type = {}' is not valid (Note: Triple frequency solution is not " "in use.).",
            freq_type,
        )

    return code_obs


# TODO: Connect needed between station and satellite position
#      Already part of Position library: posistion.distance / constant.c
def get_flight_time(dset):
    """Get flight time of GNSS signal between satellite and receiver

    Args:
        dset(Dataset):    Model data

    Return:
        numpy.ndarray:    Flight time of GNSS signal between satellite and receiver in [s]
    """
    from where.models.delay import gnss_range  # Local import to avoid cyclical import

    # Get geometric range between satellite and receiver position
    geometric_range = gnss_range.gnss_range(dset)

    return geometric_range / constant.c


def get_number_of_satellites(dset: "Dataset") -> np.ndarray:
    """Get number of satellites per epoch

    Args:
        dset:  A dataset containing the data.

    Returns:
        Number of satellites per epoch
    """
    num_satellite = np.zeros((dset.num_obs))

    for sys in dset.unique("system"):
        idx_sys = dset.filter(system=sys)
        num_satellite_epoch = np.zeros((len(dset.system[idx_sys])))

        for time in dset.unique("time"):
            idx = dset.time.datetime[idx_sys] == time.datetime
            num_satellite_epoch[idx] = len(dset.satellite[idx_sys][idx])

        num_satellite[idx_sys] = num_satellite_epoch

    return num_satellite


def obstype_to_freq(sys, obstype):
    """Get GNSS frequency based on given GNSS observation type

    Args:
        sys(str):     GNSS identifier (e.g. 'E', 'G', ...)
        obstype(str): Observation type (e.g. 'L1', 'P1', 'C1X', ...)

    Return:
        float:    GNSS frequency in [Hz]
    """
    try:
        freq = getattr(enums, "gnss_freq_" + sys)[getattr(enums, "gnss_num2freq_" + sys)["f" + obstype[1]]]
    except KeyError:
        log.fatal(f"Frequency for GNSS '{sys}' and observation type '{obstype}' is not defined.")

    return freq


def get_initial_flight_time(dset, sat_clock_corr=None, rel_clock_corr=None):
    r"""Get initial flight time of GNSS signal between satellite and receiver

    In the following it will be described, how the satellite transmission time is determined. The GNSS receiver
    registers the observation time, i.e. when the satellite signal is tracked by the receiver. In addition the
    pseudorange :math:`P_r^s` between the satellite and the receiver is observed by the GNSS receiver. The first guess
    of time of transmission :math:`t^s` can be determined if we subtract from the receiver time :math:`t_r` the time of
    flight of the GNSS signal based on the pseudorange as follows:

    .. math::
          t_0^s  = t_r - \frac{P_r^s}{c}

    with the speed of light :math:`c` and the flight time of the GNSS signal fromt the satellite to the receiver
    :math:`\frac{P_r^s}{c}`, which is determined in this function.

    The time of satellite transmission has to be corrected like:

    .. math::
        \Delta t^s = t_0^s - \Delta t_{sv} - \Delta t_r,

    with the satellite clock correction :math:`\Delta t_{sv}`:

    .. math::
         \Delta t_{sv} = a_0 + a_1 (t_0^s) + a_2 (t_0^s)^2,

    and the relativistic correction due to orbit eccentricity :math:`\Delta t_r`.

    The satellite clock correction and the relativistic eccentricity correction are applied, if this information is
    already available by the routine call.

    Args:
        dset (Dataset):                   Model data.
        sat_clock_corr (numpy.ndarray):   Satellite clock correction
        rel_clock_corr (numpy.ndarray):   Relativistic clock correction due to orbit eccentricity corrections for each
                                          observation
    Return:
       TimeDelta: Flight time of GNSS signal between satellite and receiver
    """
    # Note: It can be that the observation table 'obs' is not given. For example if different orbit solutions are
    #       compared, it is not necessary to read GNSS observation data. In this case the Dataset time entries
    #       are not corrected for time of flight determined based on pseudorange observations. Instead the given
    #       Dataset time entries are directly used.
    flight_time = np.zeros(dset.num_obs)
    if "obs" in dset.fields:
        for sys in dset.unique("system"):

            # Get code observation type defined by given observation and observation type priority list
            # Note: First element of GNSS observation type list should be used.
            obstype = dset.meta["obstypes"][sys][0]
            log.debug(
                f"Code observation '{obstype}' for GNSS '{sys}' is selected for determination of initial flight time."
            )

            idx = dset.filter(system=sys)
            flight_time[idx] = dset.obs[obstype][idx] / constant.c

    if sat_clock_corr is not None:
        flight_time += sat_clock_corr / constant.c

    if rel_clock_corr is not None:
        flight_time += rel_clock_corr / constant.c

    return TimeDelta(flight_time, fmt="seconds", scale="gps")


# TODO: already in Position via 'direction'
def get_line_of_sight(dset):
    """Get the Line of Sight vector from receiver to satellite in the ITRS.
    """
    # TODO: Other solution dset.site_pos.convert_gcrs_to_itrs(dset.site_pos.direction)
    return mathp.unit_vector(dset.sat_posvel.itrs_pos - dset.site_pos.itrs)


def get_rinex_file_version(file_key=None, file_vars=None, file_path=None):
    """ Get RINEX file version for a given file key

    Args:
        file_key:       File key defined in files.conf file (e.g. given for RINEX navigation or observation file)
        file_vars:      Variables needed to identify RINEX file based on definition in files.conf file.
        file_path (pathlib.PosixPath):  File path to broadcast orbit file.
        
    Returns:
        tuple:         with following elements

    ===============  ==================================================================================
     Elements          Description
    ===============  ==================================================================================
     version          RINEX file version
     filepath         RINEX file path
    ===============  ==================================================================================
    """
    file_vars = dict() if file_vars is None else file_vars
    if file_path is None:
        file_path = config.files.path(file_key, file_vars=file_vars)

    with config.files.open_path(file_path, mode="rt") as infile:
        try:
            version = infile.readline().split()[0]
        except IndexError:
            log.fatal(f"Could not find Rinex version in file {file_path}")

    return version, file_path


def gpssec2jd(wwww, sec):
    """
    FUNCTION: gpsSec2jd(wwww,sec)

    PURPOSE:  Conversion from GPS week and second to Julian Date (JD)

    RETURN:   (float) jd_day, jd_frac - Julian Day and fractional part

    INPUT:    (int) wwww, (float) sec - GPS week and second
    """
    SEC_OF_DAY = 86400.0
    JD_1980_01_06 = 2_444_244  # Julian date of 6-Jan-1980 + 0.5 d

    # .. Determine GPS day
    wd = np.floor((sec + 43200.0) / 3600.0 / 24.0)  # 0.5 d = 43200.0 s

    # .. Determine remainder
    fracSec = sec + 43200.0 - wd * 3600.0 * 24.0

    # .. Conversion GPS week and day to from Julian Date (JD)
    jd_day = wwww * 7.0 + wd + JD_1980_01_06
    jd_frac = fracSec / SEC_OF_DAY

    return jd_day, jd_frac


def jd2gps(jd):
    """
    FUNCTION: jd2gps(jd)

    PURPOSE:  Conversion from Julian Date (JD) to GPS week and day (started 6-Jan-1980).

    RETURN:   (int) wwww, wd, frac - GPS week, GPS day and fractional part / GPS seconds

    INPUT:    (float) jd - Julian Date
    """
    JD_1980_01_06 = 2_444_244.5  # Julian date of 6-Jan-1980
    if np.any(jd < JD_1980_01_06):
        log.fatal("Julian Day exceeds the GPS time start date of 6-Jan-1980 (JD 2444244.5).")

    # .. Conversion from Julian Date (JD) to GPS week and day
    wwww = np.floor((jd - JD_1980_01_06) / 7)
    wd = np.floor(jd - JD_1980_01_06 - wwww * 7)
    frac = jd - JD_1980_01_06 - wwww * 7 - wd
    gpssec = (frac + wd) * 86400.0

    return wwww, wd, frac, gpssec


def linear_combination(type_: str, dset: "Dataset") -> Tuple[np.ndarray]:
    """Calculate linear combination of observations for given linear combination type

    Args:
        dset:    Dataset
        type_:   Type of linear combination: 'ionosphere-free'

    Returns:
        Tuple  with following `numpy.ndarray` arrays

    ===================  ============================================================================================
     Elements             Description
    ===================  ============================================================================================
     code_obs_combined    Array with combined code observations in [m].
     phase_obs_combined   Array with combined carrier phase observations in [m].
    ===================  ============================================================================================
    """
    code_obs_combined = np.zeros(dset.num_obs)
    phase_obs_combined = np.zeros(dset.num_obs)

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)

        # Get pseudorange and carrier phase observations for the 1st and 2nd frequency
        #
        # NOTE: The GNSS observation types defined in meta variable 'obstypes' has a defined order, which is determined
        #       by the given observation types for each GNSS and the priority list.
        #
        #
        observation_code = config.tech.gnss_select_obs.obs_code.str
        if observation_code == "code":
            C1 = dset.meta["obstypes"][sys][0]  # Pseudorange observation for 1st frequency
            C2 = dset.meta["obstypes"][sys][1]  # Pseudorange observation for 2nd frequency
        elif observation_code == "phase":
            L1 = dset.meta["obstypes"][sys][0]  # Carrier phase observation for 1st frequency
            L2 = dset.meta["obstypes"][sys][1]  # Carrier phase observation for 2nd frequency
        elif observation_code == "code:phase":
            C1 = dset.meta["obstypes"][sys][0]  # Pseudorange observation for 1st frequency
            L1 = dset.meta["obstypes"][sys][1]  # Carrier phase observation for 1st frequency
            C2 = dset.meta["obstypes"][sys][2]  # Pseudorange observation for 2nd frequency
            L2 = dset.meta["obstypes"][sys][3]  # Carrier phase observation for 2nd frequency
        else:
            log.fatal(f"Linear combination determination is not defined for observation code {observation_code}.")

        if type_ == "ionosphere-free":
            if "code" in observation_code:
                f1 = getattr(enums, "gnss_freq_" + sys)["f" + C1[1]]  # Frequency of 1st band
                f2 = getattr(enums, "gnss_freq_" + sys)["f" + C2[1]]  # Frequency of 2nd band
                code_obs_combined[idx] = ionosphere_free_linear_combination(dset[C1][idx], dset[C2][idx], f1, f2)
            elif "phase" in observation_code:
                f1 = getattr(enums, "gnss_freq_" + sys)["f" + L1[1]]  # Frequency of 1st band
                f2 = getattr(enums, "gnss_freq_" + sys)["f" + L2[1]]  # Frequency of 2nd band
                phase_obs_combind[idx] = ionosphere_free_linear_combination(dset[L1][idx], dset[L2][idx], f1, f2)
        else:
            log.fatal(f"Linear combination type '{type_}' is not defined.")

    return code_obs_combined, phase_obs_combined


def ionosphere_free_linear_combination(obs1: np.ndarray, obs2: np.ndarray, f1: float, f2: float) -> np.ndarray:
    """Generate ionosphere-free linear combination

    Args:
        obs1:   1st observation array 
        obs2:   2nd observation array
        f1:     frequency of 1st observation
        f2:     frequency of 2nd observation
    """

    # Coefficient of ionospheric-free linear combination
    n = f1 ** 2 / (f1 ** 2 - f2 ** 2)
    m = -f2 ** 2 / (f1 ** 2 - f2 ** 2)

    # Generate ionospheric-free linear combination
    return n * obs1 + m * obs2


# TODO hjegei: Better solution?
def llh2xyz(lat, lon, h):
    """Conversion of geodetic (geographical) to cartesian to geodetic.

    Reference: "Geocentric Datum of Australia", Technical Manual, Version 2.4, Intergovernmental Committee on Surveying
               and Mapping, 2 December 2014, page 33

    Args:    (float) lat,lon,h  - Geodetic (geographical) coordinates latitude, east longitude in radian and ellipsoidal
                                  height in meter

    Returns:   (float) x,y,z      - Geocentric cartesian coordiantes [m]
    """

    # .. Local variables
    SEMI_MAJOR_AXIS_WGS84 = 6_378_137.0
    FLATTENING_WGS84 = 1.0 / 298.257_223_563
    a = SEMI_MAJOR_AXIS_WGS84
    f = FLATTENING_WGS84

    # .. Calculate help parameters
    e2 = (2 - f) * f  # squared eccentricity
    sin2lat = np.sin(lat) * np.sin(lat)
    v = a / np.sqrt(1 - e2 * sin2lat)

    # .. Calculate coordinates
    x = (v + h) * np.cos(lat) * np.cos(lon)
    y = (v + h) * np.cos(lat) * np.sin(lon)
    z = ((1 - e2) * v + h) * np.sin(lat)

    # .. Return geodetic coordinates in [m]
    return x, y, z


def plot_skyplot(dset):
    """Plot skyplot

    Args:
        dset
    """

    cm = plt.get_cmap("gist_rainbow")
    ax = plt.subplot(111, projection="polar")
    ax.set_prop_cycle(
        plt.cycler(
            "color", [cm(1.0 * i / len(dset.unique("satellite"))) for i in range(len(dset.unique("satellite")))]
        )
    )
    for sat in dset.unique("satellite"):
        idx = dset.filter(satellite=sat)
        azimuth = dset.site_pos.azimuth[idx]
        zenith_distance = np.rad2deg(dset.site_pos.zenith_distance[idx])
        ax.plot(azimuth, zenith_distance, ".", markersize=7, label=sat)
        ax.set_ylim(90)  # set radius of circle to the maximum elevation
        ax.set_theta_zero_location("N")  # sets 0(deg) to North
        ax.set_theta_direction(-1)  # sets plot clockwise
        ax.set_yticks(range(0, 90, 30))  # sets 3 concentric circles
        ax.set_yticklabels(map(str, range(90, 0, -30)))  # reverse labels
        ax.grid("on")
        ax.legend(fontsize=8, loc=1, bbox_to_anchor=(1.25, 1.0))
        ax.set_title("Skyplot", va="bottom")
    plt.show()
