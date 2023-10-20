"""Clean Galileo HAS orbit. Remove GNSS observations with unavailable HAS message data, and
observation exceeding the validity length of a navigation record.

Description:
------------

The remover can be used for Galileo HAS orbits to remove GNSS observations. Following steps are carried out:
    - removing of GNSS observations if HAS messages are unavailable
    - removing of GNSS observations with timestamps less than HAS receiver reception time
    - removing of GNSS observations which exceeds the validity length of a Galileo HAS record
"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import cleaners
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_clean_orbit_has(dset: "Dataset", orbit: "HasOrbit") -> None:
    """Remove GNSS observations for unavailable Galileo HAS satellite orbits or which does not fulfill requirements

    Args:keep_idx = np.ones(dset.num_obs, dtype=bool)
        dset:           A Dataset containing model data.
        orbit:          Object of HasOrbit class.

    Returns:
        numpy.ndarray:   Array containing False for observations to throw away.
    """   
    check_has_validity_length = config.tech[_SECTION].check_has_validity_length.bool

    # GNSS satellite observations are rejected from Dataset 'dset', if HAS message are not avaibale for a satellite
    _ignore_satellites(dset, orbit)
    
    # Remove GNSS observations which has no valid HAS message nearby
    _ignore_epochs_with_no_valid_has_message(dset, orbit)

    # Remove GNSS observations which exceeds the validity length of HAS messages
    keep_idx = np.ones(dset.num_obs, dtype=bool)    
    
    if check_has_validity_length:
        keep_idx = _ignore_epochs_exceeding_validity(dset, orbit)
 
    return keep_idx


def _ignore_satellites(dset: "Dataset", orbit: "Dataset") -> None:
    """Remove GNSS observations with unavailable HAS messages from Dataset

    Args:
        dset:       A Dataset containing model data, which is decimated by unavailable satellite observations.
        orbit:      Dataset containing orbit data.
    """

    # Account only for satellites missing for rundate
    idx = orbit.dset_edit.time.gps.date == orbit.dset_edit.analysis["rundate"].strftime("%Y-%m-%d")

    not_available_sat = set(dset.satellite) - set(orbit.dset_edit.satellite[idx])
    file_paths = orbit.dset_edit.meta["parser"]["file_path"]

    if not_available_sat:
        log.warn(
            f"For following satellites are no valid HAS messages available {', '.join(file_paths)}: "
            f"{', '.join(sorted(not_available_sat))}"
        )
        cleaners.apply_remover("ignore_satellite", dset, satellites=not_available_sat)
        
        
def _ignore_epochs_with_no_valid_has_message(dset: "Dataset", orbit: "Dataset") -> np.ndarray:
    """Remove GNSS observations for observation epoch with timestamps less than HAS message receiver reception time

    Args:
        dset:   A Dataset containing model data.
        orbit:      Dataset containing orbit data.
    """
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    obs_epoch_nearest_positive = True if "positive" in config.tech.has_message_nearest_to.str else False

    # Make lookup table of gps seconds timestamps of the dataset for each satellite
    gpssecs_by_sat = dict()
    for sat in dset.unique("satellite"):
        idx = orbit.dset_edit.filter(satellite=sat)
        gpssecs = orbit.dset_edit.time.gps.gps_ws.seconds[idx]
        gpssecs_by_sat.update({sat:(idx, gpssecs)})

    for i, (epoch, sat) in enumerate(zip(dset.time, dset.satellite)):

        if obs_epoch_nearest_positive:
            
            idx = gpssecs_by_sat[sat][0]
            diff = epoch.gps.gps_ws.seconds - gpssecs_by_sat[sat][1]
    
            if np.all(diff < 0):
                log.debug(f"No valid HAS message could be found for satellite {sat} and observation epoch {epoch.datetime.strftime('%Y-%d-%mT%H:%M:%S')} ")
                #log.debug(f"No valid HAS message could be found for satellite {sat} and observation epoch {epoch.datetime.strftime('%Y-%d-%mT%H:%M:%S')} "  # slow
                #          f"(nearest receiver reception time of HAS message: {min(orbit.dset_edit.time.gps.datetime[idx]).strftime('%Y-%d-%mT%H:%M:%S')},"   # slow
                #          f"{min(orbit.dset_edit.time.gps_ws.week[idx]):.0f}-{min(orbit.dset_edit.time.gps_ws.seconds[idx]):.0f})")  # old code
                keep_idx[i] = False
            
    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations based on _ignore_epochs (file key '{orbit.file_key}')")
    
    # Note: Observations have to be removed already here, otherwise further processing does not work.
    if num_removed_obs > 0:
        log.info(f"Keeping {sum(keep_idx)} of {dset.num_obs} observations.")
        dset.subset(keep_idx)
 

def _ignore_epochs_exceeding_validity(dset: "Dataset", orbit: "HasOrbit") -> np.ndarray:
    """Remove GNSS observations which exceeds the validity length of a HAS message

    Each HAS message has information about the validity length. If this validity limit is exceeded for an observation, 
    then this observation is rejected.

    Args:
        dset:   A Dataset containing model data.
        orbit:  Object of HasOrbit class.
        
    Returns:
        Array containing False for observations to throw away
    """
    has_message_idx = orbit._get_has_message_idx(dset)
    keep_idx = np.ones(dset.num_obs, dtype=bool)
    
    # Age of HAS messages
    age_of_has_message = (dset.time.gps.mjd - orbit.dset_edit.tom.gps.mjd[has_message_idx]) * Unit.day2second
    
    # Remove observations, if they exceed validity length of HAS messages
    idx = abs(age_of_has_message) > orbit.dset_edit.validity[has_message_idx]
    keep_idx[idx] = False

    ##+DEBUG
    #rejected_obs = [f"{s}, {i}, {t}, {to}, {a:.0f} > {v}\n" for s, i, t, to, a, v in zip(
    #            dset.time.gps.datetime[~keep_idx], 
    #            dset.satellite[~keep_idx], 
    #            orbit.dset_edit.iod[has_message_idx][~keep_idx],
    #            orbit.dset_edit.tom.gps.datetime[has_message_idx][~keep_idx],
    #            abs(age_of_has_message[~keep_idx]),
    #            orbit.dset_edit.validity[has_message_idx][~keep_idx]
    #)]
    #print(f"REJECT:{''.join(rejected_obs)}")
    ##-DEBUG

    num_removed_obs = dset.num_obs - np.count_nonzero(keep_idx)
    log.info(f"Removing {num_removed_obs} observations based on _ignore_epochs_exceeding_validity (file key '{orbit.file_key}')")

    #log.debug('Following entries are removed: {}\n', 'DEBUG:  \n'.join([s+t.strftime('  %Y-%m-%d %H:%M:%S (GPS)')
    #                                                              for s, t in zip(np.array(dset.satellite)[keep_idx],
    #                                                                              dset.time.gps.datetime[keep_idx])]))

    return keep_idx
