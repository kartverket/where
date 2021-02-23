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
get_earth_rotation(dset)                       -> Midgard: PosVel (see function)
get_observation(dset)                          -> Midgard: gnss
get_flight_time(dset)                          -> Midgard: PosVel (see function)
obstype_to_freq(sys, obstype)                  -> Midgard: gnss
get_initial_flight_time(dset, sat_clock_corr=None, rel_clock_corr=None)  -> Midgard: gnss
get_line_of_sight(dset)                        -> Midgard: Position library
get_rinex_file_version(file_key, file_vars)    -> Is that needed in the future?
linear_combination(dset)                       -> Midgard: gnss
llh2xyz(lat, lon, h)                           -> midgard.math.transformation.llh2trs
plot_skyplot(dset)                             -> Midgard: plot

Should we more specific in using arguments, that instead of using 'dset'? -> Maybe

"""
# Standard library imports
from typing import Any, Dict, Tuple

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


def get_observation(dset):
    """Get observations depending on given observation types

    The first element of the observation type variable `dset.meta['obstypes'][sys]` is selected as observation for
    single frequency solution. The order of the observation type variable `dset.meta['obstypes'][sys]` depends on
    the priority list given in the configuration file and the given observations.

    The ionospheric-free linear combination is applied for dual frequency solution.

    Args:
        dset:    Dataset

    Returns:
        numpy.ndarray:  Observation choosen depending on priority list and used pipeline. For dual frequency solution
                        observations are given as ionospheric-free linear combination.
    """
    # TODO: Function has to be more generalized. At the moment only SPV and SPP solutions can be handled.
    freq_type = config.tech.freq_type.str
    observation = np.zeros(dset.num_obs)

    if freq_type == "single":

        if dset.vars["pipeline"] == "gnss_vel":
            for sys in dset.unique("system"):
                idx = dset.filter(system=sys)
                obstypes = dset.meta["obstypes"][sys]

                # Note: SPV solution needs pseudo-range observations for determination of flight time between satellite
                #       and receiver. This pseudo-range observations are not used for further SPV calculation, only
                #       Doppler observations are of interest.
                obstype = None
                for type_ in obstypes:
                    if type_.startswith("D"):
                        obstype = type_
                        break

                if not obstype:
                    log.fatal(f"No Doppler observations are defined for {sys}: {', '.join(obstypes)}")

                observation[idx] = dset.obs[obstype][idx]

        else:
            for sys in dset.unique("system"):
                idx = dset.filter(system=sys)
                obstypes = dset.meta["obstypes"][sys]
                observation[idx] = dset.obs[obstypes[0]][idx]

    elif freq_type == "dual":
        linear_comb = linear_combination("ionosphere_free", dset)
        observation = linear_comb["code"]["val"]
    else:
        log.fatal(
            f"Configuration option 'freq_type = {freq_type}' is not valid (Note: Triple frequency solution is not "
            f"in use.)."
        )

    return observation


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


def obstype_to_freq(sys: str, obstype: str) -> float:
    """Get GNSS frequency based on given GNSS observation type

    Args:
        sys:     GNSS identifier (e.g. 'E', 'G', ...)
        obstype: Observation type (e.g. 'L1', 'P1', 'C1X', ...)

    Return:
        GNSS frequency in [Hz]
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


def linear_combination(type_: str, dset: "Dataset") -> Dict[str, Dict[str, Any]]:
    """Calculate linear combination of observations for given linear combination type and same observation type
    (code, phase, doppler, snr)

    Args:
        dset:    Dataset
        type_:   Type of linear combination, which can be 'geometry_free', 'ionosphere_free', 'narrow_lane' or 
                 'wide_lane'.

    Returns:
        Dictionary with observation type as key (code, phase, doppler and/or snr) and dictionary with array with 
        linear combination values as values in [m] and name of combined observations.
    """
    func = {
        "geometry_free": geometry_free_linear_combination,
        "ionosphere_free": ionosphere_free_linear_combination,
        "narrow_lane": narrowlane_linear_combination,
        "wide_lane": widelane_linear_combination,
    }

    cfg_obs_code = config.tech.gnss_select_obs.obs_code.list
    linear_comb = dict()
    for obs_code in cfg_obs_code:
        linear_comb[obs_code] = dict(val = np.zeros(dset.num_obs))

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)

        # Get observations for the 1st and 2nd frequency
        #
        # NOTE: The GNSS observation types defined in meta variable 'obstypes' has a defined order, which is determined
        #       by the given observation types for each GNSS and the priority list.
        #
        obs_num = 0
        for obs_code in cfg_obs_code:

            obs_1 = dset.meta["obstypes"][sys][obs_num]
            obs_2 = dset.meta["obstypes"][sys][obs_num + 1]
            linear_comb[obs_code].setdefault("sys_obs", dict()).update({sys: [obs_1, obs_2]})

            log.debug(
                f"Generate {type_} combination for GNSS '{sys}' and {obs_code} observations {obs_1} and {obs_2}."
            )

            if type_ == "geometry_free":
                linear_comb[obs_code]["val"][idx] = func[type_](dset.obs[obs_1][idx], dset.obs[obs_2][idx])
            else:
                f1 = getattr(enums, "gnss_freq_" + sys)["f" + obs_1[1]]  # Frequency of 1st band
                f2 = getattr(enums, "gnss_freq_" + sys)["f" + obs_2[1]]  # Frequency of 2nd band
                log.debug(
                    f"Frequencies for {type_} combination: f1 = {f1} Hz ({obs_1}), f2 = {f2} Hz ({obs_2})."
                )

                try:
                    linear_comb[obs_code]["val"][idx] = func[type_](dset.obs[obs_1][idx], dset.obs[obs_2][idx], f1, f2)
                except KeyError:
                    log.fatal(f"Linear combination 'type_' is not defined.")

            obs_num += 2

    return linear_comb


def linear_combination_cmc(dset: "Dataset") -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Calculate code multipath linear combination (CMC) based on code and phase observations

    Args:
        dset:    Dataset

    Returns:
        Tuple with code multipath linear combination for frequency 1 and 2 in [m]. Each code linear combination is
        saved in a dictionary with value of linear combination and name of combined observations.
    """
    cmc1 = dict(val = np.zeros(dset.num_obs))
    cmc2 = dict(val = np.zeros(dset.num_obs))

    for sys in dset.unique("system"):
        
        if len(dset.meta["obstypes"][sys]) < 4:
            raise ValueError(f"Dual-frequency code and phase observations are needed for code multipath linear "
                             f"combination.")
        
        idx = dset.filter(system=sys)
        
        # Get observations for the 1st and 2nd frequency
        #
        # NOTE: The GNSS observation types defined in meta variable 'obstypes' has a defined order, which is determined
        #       by the given observation types for each GNSS and the priority list.
        #
        code1 = dset.meta["obstypes"][sys][0]
        code2 = dset.meta["obstypes"][sys][1]
        phase1 = dset.meta["obstypes"][sys][2]
        phase2 = dset.meta["obstypes"][sys][3]
        cmc1.setdefault("sys_obs", dict()).update({sys: [code1, phase1, phase2]})
        cmc2.setdefault("sys_obs", dict()).update({sys: [code2, phase2, phase1]})

        f1 = getattr(enums, "gnss_freq_" + sys)["f" + code1[1]]  # Frequency of 1st band
        f2 = getattr(enums, "gnss_freq_" + sys)["f" + code2[1]]  # Frequency of 2nd band

        cmc1["val"][idx], cmc2["val"][idx] = code_multipath_linear_combination(
            dset.obs[code1][idx], dset.obs[code2][idx], dset.obs[phase1][idx], dset.obs[phase2][idx], f1, f2
        )

    return cmc1, cmc2


def code_phase_difference(dset: "Dataset") -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Calculate code-phase difference based on code and phase observations

    Args:
        dset:    Dataset

    Returns:
        Tuple with code-phase difference for frequency 1 and 2 in [m]. Each code-phase difference is saved in a 
        dictionary with value of linear combination and name of combined observations.
    """
    code_phase_1 = dict(val = np.zeros(dset.num_obs))
    code_phase_2 = dict(val = np.zeros(dset.num_obs))

    for sys in dset.unique("system"):
        
        # TODO: This is not correct for single-frequency solutions. 
        if len(dset.meta["obstypes"][sys]) < 4:
            raise ValueError(f"Dual-frequency code and phase observations are needed for code-phase difference.")
        
        idx = dset.filter(system=sys)
        
        # Get observations for the 1st and 2nd frequency
        #
        # NOTE: The GNSS observation types defined in meta variable 'obstypes' has a defined order, which is determined
        #       by the given observation types for each GNSS and the priority list.
        #
        code1 = dset.meta["obstypes"][sys][0]
        code2 = dset.meta["obstypes"][sys][1]
        phase1 = dset.meta["obstypes"][sys][2]
        phase2 = dset.meta["obstypes"][sys][3]
        code_phase_1.setdefault("sys_obs", dict()).update({sys: [code1, phase1]})
        code_phase_2.setdefault("sys_obs", dict()).update({sys: [code2, phase2]})

        code_phase_1["val"][idx] = dset.obs[code1][idx] - dset.obs[phase1][idx]
        code_phase_2["val"][idx] = dset.obs[code1][idx] - dset.obs[phase1][idx]

    return code_phase_1, code_phase_2


def linear_combination_melbourne(dset: "Dataset") -> Dict[str, Any]:
    """Calculate Melbourne-Wübbena linear combination based on code and phase observations

    Args:
        dset:    Dataset

    Returns:
        Melbourne-Wübbena linear combination
    """
    linear_comb = dict(val = np.zeros(dset.num_obs))

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)
        
        if len(dset.meta["obstypes"][sys]) < 4:
            raise ValueError("Two code and two phase observations are needed for Melbourne-Wübbena linear combination.")

        # Get observations for the 1st and 2nd frequency
        #
        # NOTE: The GNSS observation types defined in meta variable 'obstypes' has a defined order, which is determined
        #       by the given observation types for each GNSS and the priority list.
        #
        code1 = dset.meta["obstypes"][sys][0]
        code2 = dset.meta["obstypes"][sys][1]
        phase1 = dset.meta["obstypes"][sys][2]
        phase2 = dset.meta["obstypes"][sys][3]
        linear_comb.setdefault("sys_obs", dict()).update({sys: [code1, code2, phase1, phase2]})

        f1 = getattr(enums, "gnss_freq_" + sys)["f" + code1[1]]  # Frequency of 1st band
        f2 = getattr(enums, "gnss_freq_" + sys)["f" + code2[1]]  # Frequency of 2nd band

        linear_comb["val"][idx] = melbourne_linear_combination(
            dset.obs[code1][idx], dset.obs[code2][idx], dset.obs[phase1][idx], dset.obs[phase2][idx], f1, f2
        )

    return linear_comb


def code_multipath_linear_combination(
    code1: np.ndarray, code2: np.ndarray, phase1: np.ndarray, phase2: np.ndarray, f1: float, f2: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate code multipath linear combination (CMC)

    Args:
        code1:   1st code observation array 
        code2:   2nd code observation array
        phase1:  1st phase observation array 
        phase2:  2nd phase observation array
        f1:      frequency of 1st observation
        f2:      frequency of 2nd observation
        
    Returns:
        Tuple with code multipath linear combination for frequency 1 og 2
    """
    # Coefficient of linear combination
    def _get_coefficents(f1, f2):
        n = -(f1 ** 2 + f2 ** 2) / (f1 ** 2 - f2 ** 2)
        m = 2 * f2 ** 2 / (f1 ** 2 - f2 ** 2)
        return n, m

    def _get_cmc(code, phase1, phase2, f1, f2):
        n, m = _get_coefficents(f1, f2)
        return code + n * phase1 + m * phase2

    # Generate linear combination
    cmc1 = _get_cmc(code1, phase1, phase2, f1, f2)
    cmc2 = _get_cmc(code2, phase2, phase1, f2, f1)

    return cmc1, cmc2


def geometry_free_linear_combination(obs1: np.ndarray, obs2: np.ndarray) -> np.ndarray:
    """Generate geometry-free linear combination

    Args:
        obs1:   1st observation array 
        obs2:   2nd observation array
        
    Returns:
        Geometry-free linear combination
    """
    # Coefficient of linear combination
    n = 1
    m = -1

    # Generate linear combination
    return n * obs1 + m * obs2


def ionosphere_free_linear_combination(obs1: np.ndarray, obs2: np.ndarray, f1: float, f2: float) -> np.ndarray:
    """Generate ionosphere-free linear combination

    Args:
        obs1:   1st observation array 
        obs2:   2nd observation array
        f1:     frequency of 1st observation
        f2:     frequency of 2nd observation
        
    Returns:
        Ionosphere-free linear combination
    """
    # Coefficient of linear combination
    n = f1 ** 2 / (f1 ** 2 - f2 ** 2)
    m = -f2 ** 2 / (f1 ** 2 - f2 ** 2)

    # Generate linear combination
    return n * obs1 + m * obs2


def melbourne_linear_combination(
    code1: np.ndarray, code2: np.ndarray, phase1: np.ndarray, phase2: np.ndarray, f1: float, f2: float
) -> np.ndarray:
    """Generate Melbourne-Wübbena linear combination

    Args:
        code1:   1st code observation array 
        code2:   2nd code observation array
        phase1:  1st phase observation array 
        phase2:  2nd phase observation array
        f1:      frequency of 1st observation
        f2:      frequency of 2nd observation
        
    Returns:
        Melbourne-Wübbena linear combination
    """
    # Coefficient of linear combination
    k = f1 / (f1 - f2)
    l = -f2 / (f1 - f2)
    n = -f1 / (f1 + f2)
    m = -f2 / (f1 + f2)

    # Generate linear combination
    return k * phase1 + l * phase2 + n * code1 + m * code2


def widelane_linear_combination(obs1: np.ndarray, obs2: np.ndarray, f1: float, f2: float) -> np.ndarray:
    """Generate widelane linear combination

    Args:
        obs1:   1st observation array 
        obs2:   2nd observation array
        f1:     frequency of 1st observation
        f2:     frequency of 2nd observation
        
    Returns:
        Widelane linear combination
    """
    # Coefficient of linear combination
    n = f1 / (f1 - f2)
    m = -f2 / (f1 - f2)

    # Generate linear combination
    return n * obs1 + m * obs2


def narrowlane_linear_combination(obs1: np.ndarray, obs2: np.ndarray, f1: float, f2: float) -> np.ndarray:
    """Generate narrow-lane linear combination

    Args:
        obs1:   1st observation array 
        obs2:   2nd observation array
        f1:     frequency of 1st observation
        f2:     frequency of 2nd observation
        
    Returns:
        Narrow-lane linear combination
    """
    # Coefficient of linear combination
    n = f1 / (f1 + f2)
    m = f2 / (f1 + f2)

    # Generate linear combination
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
