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

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where import cleaners
from where.lib import config
from where.lib import log
from where.data.time import Time, TimeDelta

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_clean_orbit(dset: "Dataset", orbit_flag: str) -> None:
    """Remove GNSS observations with unavailable apriori satellite orbits or which does not fulfill requirements

    Args:
        dset:           A Dataset containing model data.
        orbit_flag:    Specification of which apriori orbit is used ("broadcast" or "precise")

    Returns:
        numpy.ndarray:   Array containing False for observations to throw away.
    """
    check_nav_validity_length = config.tech[_SECTION].check_nav_validity_length.bool
    ignore_unhealthy_satellite = config.tech[_SECTION].ignore_unhealthy_satellite.bool

    # GNSS observations are rejected from Dataset 'dset', if apriori satellite orbits are not given
    _ignore_satellites(dset, orbit_flag)

    # Remove unhealthy satellites
    if (orbit_flag == "broadcast") and ignore_unhealthy_satellite:
        cleaners.apply_remover("gnss_ignore_unhealthy_satellite", dset)

    # Remove GNSS observations which exceeds the validity length of broadcast ephemeris
    if (orbit_flag == "broadcast") and check_nav_validity_length:
        keep_idx = _ignore_epochs_exceeding_validity(dset)

    # Remove GNSS observations which exceeds the interpolation boundaries
    if orbit_flag == "precise":
        keep_idx = _ignore_epochs_exceeding_interpolation_boundaries(dset)

    return keep_idx


def _ignore_satellites(dset: "Dataset", orbit_flag: str) -> None:
    """Remove GNSS observations with unavailable apriori satellite orbits from Dataset

    The remover can be used for precise and broadcast ephemeris and also for the GNSS and SISRE technique/analysis.

    - GNSS: The apriori orbit is chosen via the configuration. The available satellites of the apriori orbits are
            compared against the satellites given in the GNSS observation files.

    - SISRE: Both apriori orbits, broadcast and precise, has to be checked against a set of satellites defined in
             the configuration file.

    Args:
        dset:        A Dataset containing model data, which is decimated by unavailable satellite observations.
        orbit_flag:  Specification of which apriori orbit is used ("broadcast" or "precise")
    """
    orbit = apriori.get(
        "orbit",
        apriori_orbit=orbit_flag,
        rundate=dset.analysis["rundate"],
        station=dset.vars["station"],
        system=tuple(dset.unique("system")),
        day_offset=0,  # check satellite availability only for current rundate and not for the days before/after
        # rundate in addition. TODO: This does not work for broadcast ephemeris at the moment.
    )
    not_available_sat = set(dset.satellite) - set(orbit.dset_raw.satellite)
    file_paths = orbit.dset_raw.meta["parser"]["file_path"]

    if not_available_sat:
        log.warn(
            f"The following satellites are not given in apriori {orbit_flag} orbit file {', '.join(file_paths)}: "
            f"{', '.join(sorted(not_available_sat))}"
        )
        cleaners.apply_remover("ignore_satellite", dset, satellites=not_available_sat)


def _ignore_epochs_exceeding_validity(dset: "Dataset") -> np.ndarray:
    """Remove GNSS observations which exceeds the validity length of a broadcast navigation record 

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
        dset:   A Dataset containing model data.
    """
    brdc = apriori.get(
        "orbit",
        rundate=dset.analysis["rundate"],
        station=dset.vars["station"],
        system=tuple(dset.unique("system")),
        apriori_orbit="broadcast",
    )

    brdc_block_idx = brdc._get_brdc_block_idx(dset)
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    # Loop over subset of Dataset, which includes only observations from available satellites
    for obs, (idx, time) in enumerate(zip(brdc_block_idx, dset.time.gps.mjd)):

        tk = (time - brdc.dset_edit.toe.gps.mjd[idx]) * Unit.day2second
        sys = np.array(dset.system)[obs]

        # Check validity length of navigation record
        if sys == "C":
            # TODO: :cite:`bds-sis-icd-2.1` does not define validity length of navigation record
            fit_interval = 1.0  # Assumption due to update rate of ephemeris of 1 hours
            toe_limit = fit_interval * 3600.0

        elif sys == "E":
            # Galileo navigation data record is valid for 4 hours after time of ephemeris due to Appendix C.4.4.1 in
            # Galileo-OS-SDD (2016).
            fit_interval = 4.0
            toe_limit = fit_interval * 3600.0

            # Only observation epochs after time of ephemeris should be used for Galileo. Therefore epochs with negative
            # tk has to be removed.
            if tk < 0:
                keep_idx[obs] = False
                continue

        elif sys == "G":

            # TODO: Due to :cite:`rinex3`, was the implementation of the fit interval field from the GPS navigation message
            #       an issue for the RINEX 3.02. Some implementations wrote the flag and others wrote a time interval.
            #       The RINEX 3.03 release specifies that the fit interval should be a time period for GPS and a flag
            #       for QZSS. TPP navigation files write a flag instead of a time interval for GPS, whereby 0 = 4h and
            #       1 > 4h. Should it be handled in the RINEX parser?

            # GPS navigation data record is valid for (TOE - fit_interval/2 <= epoch < TOE + fit_interval/2)
            fit_interval = 4.0 if brdc.dset_edit.fit_interval[idx] == 0.0 else brdc.dset_edit.fit_interval[idx]
            toe_limit = fit_interval * 1800.0  # toe_limit = fit_interval/2 * 3600 = fit_interval * 1800

        elif sys == "I":
            # TODO: :cite:`irnss-icd-sps` does not define validity length of navigation record
            fit_interval = 2.0  # Assumption due to update rate of ephemeris of 2 hours
            toe_limit = fit_interval * 3600.0

        elif sys == "J":
            # TODO: Due to :cite:`rinex3`, was the implementation of the fit interval field from the GPS navigation message
            #       an issue for the RINEX 3.02. Some implementations wrote the flag and others wrote a time interval.
            #       The RINEX 3.03 release specifies that the fit interval should be a time period for GPS and a flag
            #       for QZSS. TPP navigation files write a flag instead of a time interval for GPS, whereby 0 = 4h and
            #       1 > 4h. Should it be handled in the RINEX parser?

            # QZSS navigation data record is valid for (TOE - fit_interval/2 <= epoch < TOE + fit_interval/2) due to
            # section 4.1.2.7 in :cite:`is-qzss-pnt-001`
            fit_interval = 2.0 if brdc.dset_edit.fit_interval[idx] == 0.0 else brdc.dset_edit.fit_interval[idx]
            toe_limit = fit_interval * 1800.0  # toe_limit = fit_interval/2 * 3600 = fit_interval * 1800

        else:
            log.fatal(f"Broadcast ephemeris validity length interval is not defined for GNSS {sys!r}.")

        # Remove observations, if they exceed fit interval limit 'toe_limit'
        if abs(tk) > toe_limit:
            keep_idx[obs] = False

        ##+DEBUG
        #        print('DEBUG: {:6s} {:8d} {:4s} {:4.0f} {}  TOE({})  abs({:7.0f}) > {:6.0f}'.format('REJECT', obs,
        #              np.array(dset.satellite)[obs], brdc.dset_edit.iode[idx], dset.time.gps.datetime[obs],
        #              brdc.dset_edit.toe.gps.datetime[idx], tk, toe_limit))
        # else:
        #        print('DEBUG: {:6s} {:8d} {:4s} {:4.0f} {}  TOE({})  abs({:7.0f}) <={:6.0f}'.format('KEEP', obs,
        #              np.array(dset.satellite)[obs], brdc.dset_edit.iode[idx], dset.time.gps.datetime[obs],
        #              brdc.dset_edit.toe.gps.datetime[idx], tk, toe_limit))
        ##-DEBUG

    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations based on _ignore_epochs_exceeding_validity")

    # log.debug('Following entries are removed: {}\n', 'DEBUG:  \n'.join([s+t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')
    #                                                              for s, t in zip(np.array(dset.satellite)[keep_idx],
    #                                                                              dset.time.gps.datetime[keep_idx])]))

    return keep_idx


def _ignore_epochs_exceeding_interpolation_boundaries(dset: "Dataset") -> None:
    """Remove GNSS observations which exceeds the interpolation boundaries

    The interpolation boundaries are defined as +/- 'precise orbit epoch interval' related to GNSS observation epoch. If
    the boundaries are exceeded, then the GNSS observation epochs are removed. 

    Args:
        dset (Dataset):            A Dataset containing model data.
    """
    precise = apriori.get(
        "orbit", rundate=dset.analysis["rundate"], station=dset.vars["station"], apriori_orbit="precise"
    )

    epoch_interval = float(precise.dset_edit.meta[dset.analysis["rundate"].strftime("%Y-%m-%d")]["epoch_interval"])
    precise_idx = precise._get_nearest_sample_point(dset.satellite, dset.time)
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    # Check if observation epochs exceed the epoch interval of the precise orbits
    diff_time = (dset.time.gps.mjd - precise.dset_edit.time.gps.mjd[precise_idx]) * Unit.day2second
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
    keep_first_epoch, first_epoch_idx = _check_first_epoch_sample_point(dset, precise, epoch_interval)
    keep_idx[first_epoch_idx] = keep_first_epoch

    keep_last_epoch, last_epoch_idx = _check_last_epoch_sample_point(dset, precise, epoch_interval)
    keep_idx[last_epoch_idx] = keep_last_epoch

    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations based on _ignore_epochs_exceeding_interpolation_boundaries")

    return keep_idx


def _check_first_epoch_sample_point(dset: "Dataset", precise, epoch_interval):
    """Keep first observation epoch depending on existing precise orbit sample points

    Precise orbit sample points are needed to carry out interpolation of precise orbits for the first observation
    epoch. If no precise orbit sample point is available before the first satellite observation epoch, then this
    epoch will be removed for this satellite.

    Args:
        dset (Dataset):            A Dataset containing model data.
        dset_idx (numpy.ndarray):  Array containing False for observations to throw away. The array is returned by
                                   function `ignore_unavailable_orbit_satellites()`, which deletes unavailable
                                   apriori orbit satellites.
        precise (PreciseOrbit):    Precise orbit object with precise orbit information.
        epoch_interval (float):    Epoch interval of precise orbit sample points

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


def _check_last_epoch_sample_point(dset, precise, epoch_interval):
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
        dset (Dataset):            A Dataset containing model data.
        dset_idx (numpy.ndarray):  Array containing False for observations to throw away. The array is returned by
                                   function `ignore_unavailable_orbit_satellites()`, which deletes unavailable
                                   apriori orbit satellites.
        precise (PreciseOrbit):    Precise orbit object with precise orbit information.
        epoch_interval (float):    Epoch interval of precise orbit sample points

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
