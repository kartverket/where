"""Clean broadcast/precise orbit. Remove GNSS observations with unavailable satellite orbits, unhealthy satellites and
observation exceeding the validity length of a navigation record.

Description:
------------

The remover can be used for precise and broadcast ephemeris to remove GNSS observations. This can be done also for the
GNSS and SISRE analysis. Following steps are carried out:
    - removing of GNSS observations if orbits are unavailable 

and only for broadcast orbits:
    - removing of GNSS observations for satellites with unhealthy satellite status
    - removing of GNSS observations which exceeds the validity length of a broadcast navigation record

and only for precise orbits:
    - removing of GNSS observations which exceeds the interpolation boundaries
"""

# Standard library imports
from typing import Union, Tuple

# External library imports
import numpy as np

# Midgard imports
from midgard.data.time import Time
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori, cleaners
from where.lib import config, log
from where.data.time import TimeDelta

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_clean_orbit(dset: "Dataset", orbit_flag: str, orbit: Union["AprioriOrbit", None] = None) -> None:
    """Remove GNSS observations with unavailable apriori satellite orbits or which does not fulfill requirements

    Args:
        dset:          A Dataset containing model data.
        orbit_flag:    Specification of which apriori orbit is used ("broadcast" or "precise")
        orbit:         Apriori orbit object containing orbit data (either broadcast or precise)

    Returns:
        numpy.ndarray:   Array containing False for observations to throw away.
    """
    apply_has_correction = config.tech.get("apply_has_correction", default=False).bool

    # Get orbit data, if not given already
    if not orbit:
        orbit = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            system=tuple(dset.unique("system")), 
            station=dset.vars["station"],
            #day_offset=0,  # check satellite availability only for current rundate and not for the days before/after
            apriori_orbit=orbit_flag,
        )

    # GNSS observations are rejected from Dataset 'dset', if apriori satellite orbits are not given
    _ignore_satellites(dset, orbit, orbit_flag)

    # Add 'navigation_idx' field to dataset and remove epochs with no valid navigation messages
    if (orbit_flag == "broadcast"):
        keep_idx = _get_navigation_records_and_ignore_epochs(dset, orbit)  

    # GNSS observation epochs are rejected with no corresponding IOD between HAS messages and broadcast navigation 
    # messages and observation with timestamps less than HAS receiver reception time
    if (orbit_flag == "broadcast") and apply_has_correction:
        keep_idx = _ignore_epochs_has(dset, orbit)

    # Remove GNSS observations which exceeds the interpolation boundaries
    if orbit_flag == "precise":
        keep_idx = _ignore_epochs_exceeding_interpolation_boundaries(dset, orbit)


    return keep_idx

#
# AUXILIARY FUNCTIONS
#

def _check_first_epoch_sample_point(
            dset: "Dataset", 
            precise: "PreciseOrbit", 
            epoch_interval: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Keep first observation epoch depending on existing precise orbit sample points

    Precise orbit sample points are needed to carry out interpolation of precise orbits for the first observation
    epoch. If no precise orbit sample point is available before the first satellite observation epoch, then this
    epoch will be removed for this satellite.

    Args:
        dset:            A Dataset containing model data.
        precise:         Precise orbit object with precise orbit information.
        epoch_interval:  Epoch interval of precise orbit sample points

    Returns:
        tuple: Tuple with array containing False for first observations to throw away and indices indicating first
               observation epoch.
    """

    # Get indices for first observation epochs
    first_idx = 0
    first_epoch_idx = np.ones(dset.num_obs, dtype=bool)
    first_epoch_idx = dset.time.gps.mjd == dset.time.gps.mjd[first_idx]

    # Get set with satellite and time entries for getting corresponding precise orbit sample points
    satellites = dset.satellite[first_epoch_idx]
    time = Time(val=dset.time.gps.datetime[first_epoch_idx], fmt="datetime", scale=dset.time.scale) - TimeDelta(
        epoch_interval, fmt="seconds", scale=dset.time.scale
    )
    precise_idx = precise._get_nearest_sample_point(satellites, time)

    # Keep observations epochs, where a precise orbit sample point exists before the first observation epoch
    diff_time = (dset.time.gps.mjd[first_epoch_idx] - precise.dset_edit.time.gps.mjd[precise_idx]) * Unit.day2second
    keep_idx = np.logical_and(diff_time < (epoch_interval + 1), diff_time > 0)

    removed_entries = "DEBUG: ".join(
        [
            f"{s} {t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')}, dt = {dt:8.2f} s (0 < dt < {epoch_interval + 1})\n"
            for s, t, dt in zip(
                satellites[np.logical_not(keep_idx)],
                dset.time.gps.datetime[first_epoch_idx][np.logical_not(keep_idx)],
                diff_time[np.logical_not(keep_idx)],
            )
        ]
    )
    log.debug(f"Following first epoch entries are removed: \n{removed_entries}")

    return keep_idx, first_epoch_idx


def _check_last_epoch_sample_point(
            dset: "Dataset", 
            precise: "PreciseOrbit", 
            epoch_interval: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Keep last observation epoch depending on existing precise orbit sample points

    Precise orbit sample points are needed to carry out interpolation of precise orbits for the last observation
    epochs. If no precise orbit sample point is available after the last satellite observation epochs, then this
    epochs will be removed for this satellite.

    The time difference between the last observation epochs and the next precise orbit sample point is determined. 'Last
    observation epoch' + 'sampling rate' is chosen as reference time for the selection of the nearest orbit sample
    point, which corresponds normally to 0:00 GPS time. If the time difference exceeds the following intervall, then the
    observation epochs are rejected:
                       -(precise orbit epoch interval + 1) < time difference < 0

    Args:
        dset:            A Dataset containing model data.
        precise:         Precise orbit object with precise orbit information.
        epoch_interval:  Epoch interval of precise orbit sample points

    Returns:
        tuple: Tuple with array containing False for last observations to throw away and indices indicating last
               observation epoch.
    """
    sampling_rate = config.tech.sampling_rate.float

    # Get indices for last observation epochs
    last_idx = -1
    last_epoch_idx = np.ones(dset.num_obs, dtype=bool)
    last_epoch_idx = (
        dset.time.gps.mjd >= dset.time.gps.mjd[last_idx] - (epoch_interval - sampling_rate) * Unit.second2day
    )

    # Get set with satellite and time entries for getting corresponding precise orbit sample points
    # Note: Sample point reference time is 'last observation epoch' + 'sampling rate', which corresponds normally to
    #       0:00 GPS time.
    satellites = dset.satellite[last_epoch_idx]
    time = Time(val=dset.time.gps.datetime[last_idx], fmt="datetime", scale=dset.time.scale) + TimeDelta(
        sampling_rate, fmt="seconds", scale=dset.time.scale
    )

    precise_idx = precise._get_nearest_sample_point(satellites, time)

    # Keep observations epochs, where a precise orbit sample point exists after the last observation epoch
    diff_time = (dset.time.gps.mjd[last_epoch_idx] - precise.dset_edit.time.gps.mjd[precise_idx]) * Unit.day2second
    keep_idx = np.logical_and(diff_time > -(epoch_interval + 1), diff_time < 0)

    removed_entries = "DEBUG: ".join(
        [
            f"{s} {t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')}, dt = {dt:8.2f} s ({-(epoch_interval + 1)} < dt < 0)\n"
            for s, t, dt in zip(
                satellites[np.logical_not(keep_idx)],
                dset.time.gps.datetime[last_epoch_idx][np.logical_not(keep_idx)],
                diff_time[np.logical_not(keep_idx)],
            )
        ]
    )
    log.debug(f"Following last epoch entries are removed: \n{removed_entries}")

    return keep_idx, last_epoch_idx


def _get_navigation_records_and_ignore_epochs(dset: "Dataset", orbit: "AprioriOrbit") -> np.ndarray:
    """Get corresponding navigation record in relation to observation epoch and remove GNSS observations after defined
    criteria 

    'navigation_idx' field is added to dataset, which is the index to broadcast navigation record in use for each 
    observation epoch.

    Following criteria are used to reject GNSS observations in relation to broadcast navigation records:
        1. GNSS observations for which no valid broadcast navigation message are available
        2. GNSS observations which exceeds the validity length of a broadcast navigation record 
        3. GNSS observations given for unhealthy satellites

    Args:
        dset:   A Dataset containing model data.
        orbit:  Apriori orbit object containing orbit data (either broadcast or precise)
        
    Returns:
        Array containing False for observations to throw away 
    """
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    dset_brdc_idx = np.zeros(dset.num_obs, dtype=int)

    check_nav_validity_length = config.tech[_SECTION].check_nav_validity_length.bool
    ignore_unhealthy_satellite = config.tech[_SECTION].ignore_unhealthy_satellite.bool

    # Get configuration option
    brdc_block_nearest_to_options = [
        "toc",
        "toc:positive",
        "toe",
        "toe:positive",
        "transmission_time",
        "transmission_time:positive",
    ]

    brdc_block_nearest_to = config.tech.get("brdc_block_nearest_to", default="toe:positive").str.rsplit(":", 1)
    if ":".join(brdc_block_nearest_to) not in brdc_block_nearest_to_options:
        log.fatal(
            f"Unknown value {':'.join(brdc_block_nearest_to)!r} for configuration option 'brdc_block_nearest_to'. "
            f"The following values can be selected: {', '.join(brdc_block_nearest_to_options)}"
        )

    time_key = brdc_block_nearest_to[0]
    positive = True if "positive" in brdc_block_nearest_to else False

    # Loop over all observation epochs
    for obs, (satellite, time) in enumerate(zip(dset.satellite, dset.time.gps.mjd)):
        
        idx_sat = orbit.dset_edit.filter(satellite=satellite)
        idx_obs = orbit.get_nearest_idx(
                            idx = idx_sat, 
                            obs_epoch = time, 
                            time_key = time_key, 
                            positive = positive, 
                            satellite = satellite, 
                            log_level = "debug",
        )

        # Skip epochs for which no broadcast ephemeris are available
        if idx_obs is None:
            keep_idx[obs] = False
            continue  

        idx_obs = idx_sat.nonzero()[0][idx_obs]
                     
        if check_nav_validity_length:    
            status = _is_epoch_in_navigation_validity_interval(
                                time, 
                                orbit.dset_edit.toe.gps.mjd[idx_obs], 
                                orbit.dset_edit.fit_interval[idx_obs], 
                                orbit.dset_edit.iode[idx_obs],
                                satellite,
            )
            if status == False:
                keep_idx[obs] = False
                continue
        
        dset_brdc_idx[obs] = idx_obs

    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations exceeding validity length")
    log.info(f"Keeping {sum(keep_idx)} of {dset.num_obs} observations")
    # log.debug('Following entries are removed: {}\n', 'DEBUG:  \n'.join([s+t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')
    #                                                              for s, t in zip(np.array(dset.satellite)[keep_idx],
    #                                                                              dset.time.gps.datetime[keep_idx])]))

    dset.add_float("navigation_idx", val=dset_brdc_idx)
    dset.subset(keep_idx)
    
    # Remove unhealthy satellites
    if ignore_unhealthy_satellite:
        remove_idx = orbit.signal_health_status(dset) > 0

        unhealthy_satellites = sorted(set(dset.satellite[remove_idx]))
        log.info(f"Discarding observations for unhealthy satellites: {', '.join(unhealthy_satellites)}")
    else:
        remove_idx = np.zeros(dset.num_obs, dtype=bool)

    return ~remove_idx


def _get_time_of_ephemeris_limit(fit_interval: float, system: str) -> float:
    """ Get time of ephemeris limit

    How long a broadcast ephemeris block is valid depends on the GNSS:

        - BeiDou:   Not defined. Fit interval of 1 hour is used. This is an assumption due to update rate of ephemeris 
                    of 1 hours.
        - Galileo:  See appendix C.4.4.1 in :cite:`galileo-os-sdd`.
        - GPS:      Indicates the curve-fit interval used by the GPS Control Segment in determining the ephemeris
                    parameters, which is given in HOURS (see section 6.6 in :cite:`rinex2`).
        - IRNSS:    Not defined. Fit interval of 2 hour is used. This is an assumption due to update rate of ephemeris
                    of 2 hours.
        - QZSS:     Fit interval is given as flag (see section 4.1.2.7 in :cite:`is-qzss-pnt-001`):
                         0 - 2 hours
                         1 - more than 2 hours
                         blank - not known

    Args:
        fit_interval: Validity interval limit of navigation message
        system:       GNSS system

    Returns:
        Time of ephemeris limit
    """

    # Check validity length of navigation record
    if system == "C":
        # TODO: :cite:`bds-sis-icd-2.1` does not define validity length of navigation record
        fit_interval_def = 1.0  # Assumption due to update rate of ephemeris of 1 hours
        toe_limit = fit_interval_def * 3600.0

    elif system == "E":
        # Galileo navigation data record is valid for 4 hours after time of ephemeris due to Appendix C.4.4.1 in
        # Galileo-OS-SDD (2016).
        fit_interval_def = 4.0
        toe_limit = fit_interval_def * 3600.0

    elif system == "G":

        # TODO: Due to :cite:`rinex3`, was the implementation of the fit interval field from the GPS navigation message
        #       an issue for the RINEX 3.02. Some implementations wrote the flag and others wrote a time interval.
        #       The RINEX 3.03 release specifies that the fit interval should be a time period for GPS and a flag
        #       for QZSS. TPP navigation files write a flag instead of a time interval for GPS, whereby 0 = 4h and
        #       1 > 4h. Should it be handled in the RINEX parser?

        # GPS navigation data record is valid for (TOE - fit_interval/2 <= epoch < TOE + fit_interval/2)
        fit_interval = 4.0 if fit_interval == 0.0 else fit_interval
        toe_limit = fit_interval * 1800.0  # toe_limit = fit_interval/2 * 3600 = fit_interval * 1800

    elif system == "I":
        # TODO: :cite:`irnss-icd-sps` does not define validity length of navigation record
        fit_interval = 2.0  # Assumption due to update rate of ephemeris of 2 hours
        toe_limit = fit_interval * 3600.0

    elif system == "J":
        # TODO: Due to :cite:`rinex3`, was the implementation of the fit interval field from the GPS navigation message
        #       an issue for the RINEX 3.02. Some implementations wrote the flag and others wrote a time interval.
        #       The RINEX 3.03 release specifies that the fit interval should be a time period for GPS and a flag
        #       for QZSS. TPP navigation files write a flag instead of a time interval for GPS, whereby 0 = 4h and
        #       1 > 4h. Should it be handled in the RINEX parser?

        # QZSS navigation data record is valid for (TOE - fit_interval/2 <= epoch < TOE + fit_interval/2) due to
        # section 4.1.2.7 in :cite:`is-qzss-pnt-001`
        fit_interval = 2.0 if fit_interval == 0.0 else fit_interval
        toe_limit = fit_interval * 1800.0  # toe_limit = fit_interval/2 * 3600 = fit_interval * 1800

    else:
        log.fatal(f"Broadcast ephemeris validity length interval is not defined for GNSS {system!r}.")

    return toe_limit


def _ignore_epochs_exceeding_interpolation_boundaries(dset: "Dataset", orbit: "AprioriOrbit") -> np.ndarray:
    """Remove GNSS observations which exceeds the interpolation boundaries

    The interpolation boundaries are defined as +/- 'precise orbit epoch interval' related to GNSS observation epoch. If
    the boundaries are exceeded, then the GNSS observation epochs are removed. 

    Args:
        dset:   A Dataset containing model data.
        orbit:  Apriori orbit object containing orbit data (either broadcast or precise)
        
    Returns:
        Array containing False for observations to throw away 
    """
    epoch_interval = float(orbit.dset_edit.meta[dset.analysis["rundate"].strftime("%Y-%m-%d")]["epoch_interval"])
    precise_idx = orbit._get_nearest_sample_point(dset.satellite, dset.time)
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    # Check if observation epochs exceed the epoch interval of the precise orbits
    diff_time = (dset.time.gps.mjd - orbit.dset_edit.time.gps.mjd[precise_idx]) * Unit.day2second
    keep_idx = abs(diff_time) <= epoch_interval

    removed_entries = "DEBUG: ".join(
        [
            f"{s} {t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')}, diff_time: abs({dt:8.2f} s) > {epoch_interval} s\n"
            for s, t, dt in zip(
                np.array(dset.satellite[np.logical_not(keep_idx)]),
                dset.time.gps.datetime[np.logical_not(keep_idx)],
                diff_time[np.logical_not(keep_idx)],
            )
        ]
    )
    log.debug(f"Following entries are removed: \n{removed_entries}")

    # Check if first and last observation epochs exceed the epoch interval of the precise orbits
    keep_first_epoch, first_epoch_idx = _check_first_epoch_sample_point(dset, orbit, epoch_interval)
    keep_idx[first_epoch_idx] = keep_first_epoch

    keep_last_epoch, last_epoch_idx = _check_last_epoch_sample_point(dset, orbit, epoch_interval)
    keep_idx[last_epoch_idx] = keep_last_epoch

    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations based on _ignore_epochs_exceeding_interpolation_boundaries")

    # Remove all observation of a satellite, if only a single observation epoch is left (no interpolation possible)
    for sat in dset.unique("satellite"):
        sat_idx = dset.filter(satellite=sat)
        idx = dset.satellite[keep_idx] == sat
        if dset.satellite[keep_idx][idx].size == 1:
            log.warn(f"All observations of satellite {sat} are removed, because only a single observation epoch is "
                     f"left for satellite {sat} after orbit data cleaning (no interpolation possible).")
            keep_idx[sat_idx] = False
            continue

        data_period = (np.max(dset.time.gps.mjd[keep_idx][idx]) - np.min(dset.time.gps.mjd[keep_idx][idx])) * Unit.day2second
        if data_period <= epoch_interval:
            log.warn(f"All observations of satellite {sat} are removed, because the observation period of satellite "
                     f"{sat} ({data_period} s) is below orbit data interval ({epoch_interval} s) and therefore is the "
                     f" interpolation not possible.")
            keep_idx[sat_idx] = False        

    return keep_idx


def _ignore_epochs_has(dset: "Dataset", orbit: "AprioriOrbit") -> np.ndarray:
    """Remove GNSS observations for which no corresponding IOD can be found for HAS and broadcast navigation message
    and observation epoch with timestamps less than HAS message receiver reception time

    Args:
        dset:   A Dataset containing model data.
        orbit:  Apriori orbit object containing orbit data (either broadcast or precise)
        
    Returns:
        Array containing False for observations to throw away 
    """    
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    obs_epoch_nearest_positive = True if "positive" in config.tech.has_message_nearest_to.str else False
    
    for epoch, sat, iode in zip(dset.time, dset.satellite, dset.has_gnssiod_orb):

        idx = orbit.dset_edit.filter(satellite=sat, iode=iode)
        
        if np.any(idx) == False:
            idx_dset = dset.filter(satellite=sat, has_gnssiod_orb=iode)
            keep_idx[idx_dset] = False
            log.debug(f"No valid broadcast navigation message could be found for satellite {sat} and HAS message IOD "
                      f"{iode}. Removing {sum(~keep_idx[idx_dset])} observations.")
            continue
        
        if obs_epoch_nearest_positive:
           
            diff = epoch.gps.gps_ws.seconds - orbit.dset_edit.time.gps.gps_ws.seconds[idx]
    
            if np.all(diff < 0):
                log.debug(f"No valid broadcast navigation message could be found for satellite {sat}, HAS message IOD " 
                         f"{iode} and observation epoch {epoch.datetime.strftime('%Y-%d-%mT%H:%M:%S')} (nearest " 
                         f"receiver reception time of HAS message: " 
                         f"{min(orbit.dset_edit.time.gps.datetime[idx]).strftime('%Y-%d-%mT%H:%M:%S')}, "
                         f"{min(orbit.dset_edit.time.gps_ws.week[idx]):.0f}-"
                         f"{min(orbit.dset_edit.time.gps_ws.seconds[idx]):.0f})")
                idx_dset = np.logical_and(dset.filter(satellite=sat, has_gnssiod_orb=iode), epoch.gps.gps_ws.seconds == dset.time.gps.gps_ws.seconds)
                keep_idx[idx_dset] = False
            
    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations based on _ignore_epochs_has")
    
    # Note: Observations have to be removed already here, otherwise further processing does not work.
    if num_removed_obs > 0:
        log.info(f"Keeping {sum(keep_idx)} of {dset.num_obs} observations.")
        dset.subset(keep_idx)
        
    return np.ones(dset.num_obs, dtype=bool)
   

def _ignore_satellites(dset: "Dataset", orbit: "AprioriOrbit", orbit_flag: str) -> None:
    """Remove GNSS observations with unavailable apriori satellite orbits from Dataset

    The remover can be used for precise and broadcast ephemeris and also for the GNSS and SISRE technique/analysis.

    - GNSS: The apriori orbit is chosen via the configuration. The available satellites of the apriori orbits are
            compared against the satellites given in the GNSS observation files.

    - SISRE: Both apriori orbits, broadcast and precise, has to be checked against a set of satellites defined in
             the configuration file.

    Args:
        dset:        A Dataset containing model data, which is decimated by unavailable satellite observations.
        orbit:       Apriori orbit object containing orbit data (either broadcast or precise)
        orbit_flag:  Specification of which apriori orbit is used ("broadcast" or "precise")
    """
    not_available_sat = set(dset.satellite) - set(orbit.dset_edit.satellite)
    file_paths = orbit.dset_edit.meta["parser"]["file_path"]

    if not_available_sat:
        log.warn(
            f"The following satellites are not given in apriori {orbit_flag} orbit file {', '.join(file_paths)}: "
            f"{', '.join(sorted(not_available_sat))}"
        )
        cleaners.apply_remover("ignore_satellite", dset, satellites=not_available_sat)


def _is_epoch_in_navigation_validity_interval(
                    time: float, 
                    toe: float, 
                    fit_interval: float, 
                    iode: float, 
                    satellite: str, 
) -> bool:
    """Check if GNSS observation epoch is in validity length of a broadcast navigation record 

    Args:
        time:         Observation epoch in Modified Julian Day
        toe:          Time of ephemeris (reference ephemeris epoch) in Modified Julian Day
        fit_interval: Validity interval limit of navigation message
        iode:         Ephemeris issue of data (meaning depending on GNSS)
        satellite:    Satellite name

    Returns:
        True if GNSS observation epoch is in validity length of broadcast navigation record, otherwise False
    """
    brdc_block_nearest_to = config.tech.brdc_block_nearest_to.str

    tk = (time - toe) * Unit.day2second
    system = satellite[0:1]

    # Note: Epochs with negative tk has to be removed, if the difference between the observation epoch and TOE, TOC  
    #       or transmission time has to be positive. Only observation epochs after time of ephemeris should be used 
    #       for Galileo.
    if (system == "E" and tk < 0) or ("positive" in brdc_block_nearest_to and tk < 0):
        return False

    # Remove observations, if they exceed fit interval limit 'toe_limit'
    toe_limit = _get_time_of_ephemeris_limit(fit_interval, system)
    if abs(tk) > toe_limit:
        log.debug('{:6s} {:4s} {:4.0f} {}  TOE({})  abs({:7.0f}) > {:6.0f}'.format(
                            'REJECT', 
                            satellite, 
                            iode, 
                            Time(time, scale='gps', fmt='mjd').datetime,
                            Time(toe, scale='gps', fmt='mjd').datetime, 
                            tk, 
                            toe_limit,
            )
        )
        return False

    log.debug('{:6s} {:4s} {:4.0f} {}  TOE({})  abs({:7.0f}) <={:6.0f}'.format(
                        'KEEP', 
                         satellite, 
                         iode,
                         Time(time, scale='gps', fmt='mjd').datetime,
                         Time(toe, scale='gps', fmt='mjd').datetime, 
                         tk, 
                         toe_limit,
        )
    )

    return True



