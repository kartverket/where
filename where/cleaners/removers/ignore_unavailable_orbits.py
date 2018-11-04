"""Remove GNSS observations with unavailable apriori satellite orbits

Description:
------------

The remover can be used for precise and broadcast ephemeris to remove GNSS observations, where satellite orbits are
unavailable. This can be done also for the GNSS and SISRE technique/analysis. In addition can the remover remove
GNSS observations which exceeds the broadcast ephemeris fit interval.

"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.lib import plugins
from where.lib.unit import unit

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def ignore_unavailable_orbits(dset):
    """Remove GNSS observations with unavailable apriori satellite orbits

    Args:
        dset (Dataset):  A Dataset containing model data.

    Returns:
        numpy.ndarray:   Array containing False for observations to throw away.
    """
    use_remover_flag = config.tech[_SECTION].orbit.str
    orb_flag = config.tech.apriori_orbit.list

    if use_remover_flag not in ["sat", "sat_fit"]:
        return np.ones(dset.num_obs, dtype=bool)

    # Remove GNSS observations with unavailable apriori satellite orbits
    keep_idx = ignore_unavailable_apriori_orbit_satellites(dset, orb_flag)

    # Remove GNSS observations which exceeds the broadcast ephemeris fit interval
    if ("broadcast" in orb_flag) and (use_remover_flag == "sat_fit"):
        keep_idx = ignore_epochs_exceeding_brdc_fit_interval(dset, keep_idx)

    return keep_idx


def ignore_unavailable_apriori_orbit_satellites(dset, orb_flag):
    """Remove GNSS observations with unavailable apriori satellite orbits

    The remover can be used for precise and broadcast ephemeris and also for the GNSS and SISRE technique/analysis.

    - GNSS: The apriori orbit is chosen via the configuration. The available satellites of the apriori orbits are
            compared against the satellites given in the GNSS observation files.

    - SISRE: Both apriori orbits, broadcast and precise, has to be checked against a set of satellites defined in
             the configuration file.

    Args:
        dset (Dataset):   A Dataset containing model data.
        orb_flag (list):  List with used orbit types (e.g. 'broadcast' and/or 'precise' orbit)

    Returns:
        numpy.ndarray:   Array containing False for observations to throw away
    """
    remove_idx = np.zeros(dset.num_obs, dtype=bool)

    for orb in orb_flag:
        orbit = apriori.get(
            "orbit",
            apriori_orbit=orb,
            rundate=dset.rundate,
            time=dset.time,
            satellite=tuple(dset.satellite),
            system=tuple(dset.system),
            station=dset.vars["station"],
        )
        not_available_sat = set(dset.satellite) - set(orbit.dset_raw.satellite)
        file_paths = orbit.dset_raw.meta["parser"]["file_path"]

        if not_available_sat:
            log.warn(
                "Following satellites are not given in apriori {} orbit file {}: {}",
                orb,
                ", ".join(file_paths),
                ", ".join(sorted(not_available_sat)),
            )

            for satellite in not_available_sat:
                remove_idx = np.logical_or(remove_idx, dset.filter(satellite=satellite))

    num_removed_obs = str(np.count_nonzero(remove_idx))
    log.info("Removing {} observations based on ignore_unavailable_apriori_orbit_satellites", num_removed_obs)

    return np.logical_not(remove_idx)


def ignore_epochs_exceeding_brdc_fit_interval(dset, dset_idx):
    """Remove GNSS observations which exceeds the broadcast ephemeris fit interval

    How long a broadcast ephemeris block is valid depends on the GNSS:

        - GPS: Indicates the curve-fit interval used by the GPS Control Segment in determining the ephemeris
               parameters, which is given in HOURS (see section 6.6 in :cite:`rinex2`).
        - QZSS: Fit interval is given as flag (see section 4.1.2.7 in :cite:`is-qzss-pnt-001`):
                         0 - 2 hours
                         1 - more than 2 hours
                         blank - not known
        - Galileo: See appendix C.4.4.1 in :cite:`galileo-os-sdd`.

    It is necessary that the Dataset does not include unavailable apriori orbit satellites. Otherwise Where processing
    stops by getting the broadcast block indices. Therefore the argument `dset_idx` is needed, which includes
    information about unavailable apriori orbit satellites. With `dset_idx` indices the Dataset `dset` can be reduced
    to satellites, which are available in the broadcast ephemeris file.

    Args:
        dset (Dataset):            A Dataset containing model data.
        dset_idx (numpy.ndarray):  Array containing False for observations to throw away. The array is returned by
                                   function `ignore_unavailable_apriori_orbit_satellites()`, which deletes unavailable
                                   apriori orbit satellites.

    Returns:
        numpy.ndarray: Array containing False for observations to throw away
    """
    brdc = apriori.get(
        "orbit",
        rundate=dset.rundate,
        time=dset.time[dset_idx],
        satellite=tuple(np.array(dset.satellite)[dset_idx]),
        system=tuple(np.array(dset.system)[dset_idx]),
        station=dset.vars["station"],
        apriori_orbit="broadcast",
    )

    brdc_block_idx = brdc._get_brdc_block_idx()
    keep_idx = dset_idx.copy()
    subset_idx = np.ones(len(keep_idx[dset_idx]), dtype=bool)

    # Loop over subset of Dataset, which includes only observations from available satellites
    for obs, (idx, time) in enumerate(zip(brdc_block_idx, dset.time.gps.mjd[dset_idx])):
        sys = np.array(dset.system)[dset_idx][obs]
        tk = (time - brdc.dset_edit.toe.gps.mjd[idx]) * unit.day2second
        if sys == "G":

            # TODO: Due to :cite:`rinex3`, was the implementation of the fit interval field from the GPS navigation message
            #       an issue for the RINEX 3.02. Some implementations wrote the flag and others wrote a time interval.
            #       The RINEX 3.03 release specifies that the fit interval should be a time period for GPS and a flag
            #       for QZSS. TPP navigation files write a flag instead of a time interval for GPS, whereby 0 = 4h and
            #       1 > 4h. Should it be handled in the RINEX parser?
            fit_interval = 4.0 if brdc.dset_edit.fit_interval[idx] == 0.0 else brdc.dset_edit.fit_interval[idx]
            toe_limit = fit_interval * 1800.0  # toe_limit = fit_interval/2 * 3600 = fit_interval * 1800
        elif sys == "E":
            # Galileo navigation data record is valid for 4 hours after time of ephemeris due to Appendix C.4.4.1 in
            # Galileo-OS-SDD (2016).
            fit_interval = 4.0
            toe_limit = fit_interval * 3600.0

        # Remove observations, if they exceed fit interval limit 'toe_limit'
        if abs(tk) > toe_limit:
            subset_idx[obs] = False
            # TODO: This does not work keep_idx[dset_idx][obs] is still 'True'. Why?
            # keep_idx[dset_idx][obs] = False

        ##+DEBUG
        #    print('DEBUG: {:6s} {:8d} {:4s} {}  TOE({})  abs({:7.0f}) > {:6.0f}'.format('REJECT', obs,
        #          np.array(dset.satellite)[dset_idx][obs], dset.time.gps.datetime[dset_idx][obs],
        #          brdc.dset_edit.toe.gps.datetime[idx], tk, toe_limit))
        # else:
        #    print('DEBUG: {:6s} {:8d} {:4s} {}  TOE({})  abs({:7.0f}) <={:6.0f}'.format('KEEP', obs,
        #          np.array(dset.satellite)[dset_idx][obs], dset.time.gps.datetime[dset_idx][obs],
        #          brdc.dset_edit.toe.gps.datetime[idx], tk, toe_limit))
        ##-DEBUG

    # Update with keep_idx with removed observations
    for idx, value in zip(dset_idx.nonzero()[0], subset_idx):
        keep_idx[idx] = value

    num_removed_obs = np.count_nonzero(dset_idx) - np.count_nonzero(keep_idx)
    log.info("Removing {} observations based on ignore_epochs_exceeding_brdc_fit_interval", num_removed_obs)

    # log.debug('Following entries are removed: {}\n', 'DEBUG:  \n'.join([s+t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')
    #                                                              for s, t in zip(np.array(dset.satellite)[keep_idx],
    #                                                                              dset.time.gps.datetime[keep_idx])]))

    return keep_idx
