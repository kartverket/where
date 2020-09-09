"""Select GNSS observations used in Where processing

Description:
------------

    GNSS observations are selected after different criteria:

        1. Remove observations from observation types and GNSS, which are not defined in configuration file.
        2. Use of both code and phase observation or only code observation.
        3. Selection of observation types after a GNSS dependent priority list.
        4. Selection depending on using single-, dual- or triple-frequencies.

"""
# Standard library imports
import re
from typing import Set

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


OBS_CODE_DEF = dict(code="C", doppler="D", phase="L", snr="S")

FREQ_NUMBER_DEF = dict(single=["1"], dual=["1", "2"], triple=["1", "2", "3"])


# TODO MURKS: The editor with highest number is processed at the end.
@plugins.register_ordered(-100)
def gnss_select_obs(dset: "Dataset") -> np.ndarray:
    """Select GNSS observations used in Where processing

    Args:
        dset (where.data.dataset.Dataset):  A Dataset containing model data.

    Returns:
        Array containing False for observations to throw away
    """
    remove_obstypes = set()
    keep_idx = np.full(dset.num_obs, True, dtype=bool)
    obstypes_all = dset.obs.fields

    cfg_obs_code = config.tech[_SECTION].obs_code.list
    cfg_obstypes = config.tech[_SECTION].obs_types.list
    cfg_systems = config.tech.systems.list

    # Remove GNSS, which are not defined in configuration file
    for sys in list(dset.meta["obstypes"]):
        if sys not in cfg_systems:
            del dset.meta["obstypes"][sys]

    for obs, sys in enumerate(dset.system):
        if sys not in cfg_systems:
            keep_idx[obs] = False

    if not np.any(keep_idx):
        log.fatal(f"No observations available for selected system(s): {' '.join(cfg_systems)}.")

    # Remove observation types, which are not given in configuration file. If no observation types are defined in
    # configuration file keep all observation types.
    if cfg_obstypes:
        for type_ in cfg_obstypes:
            if type_ not in dset.obs.fields:
                log.warn(f"Selected observation type {type_} is not included in GNSS observation data.")
        log.debug(
            f"Remove undefined observation types in configuration file: {' '.join(set(obstypes_all) - set(cfg_obstypes))}."
        )
        remove_obstypes = set(obstypes_all) - set(cfg_obstypes)

    # Remove undefined observation codes related to given configuration
    keep_obs_code = list()
    for obs_code in sorted(cfg_obs_code):
        if obs_code not in OBS_CODE_DEF:
            log.fatal(f"Observation code '{obs_code}' is not valid in option 'obs_code='.")
        keep_obs_code.append(OBS_CODE_DEF[obs_code])

    log.debug(f"Remove undefined observation codes: {' '.join(set(OBS_CODE_DEF.values()) - set(keep_obs_code))}.")
    remove_obs_code = set(OBS_CODE_DEF.values()) - set(keep_obs_code)
    remove_obs_pattern = f"^{'|^'.join(remove_obs_code)}"

    for type_ in obstypes_all:
        search_obj = re.search(remove_obs_pattern, type_)
        if search_obj is not None:
            remove_obstypes.add(search_obj.string)

    # Select observations based on priority list
    #   -> 1st step remove already unused observation types from Dataset to determine the basis for the priority list
    #      selection

    # Note: The order of the selected observations is important for selection of GNSS code observation type to
    #       determine satellite transmission time.
    _remove_obstype_from_dset(dset, remove_obstypes)
    selected_obstypes, add_remove_obstypes = _select_observations(obstypes_all, dset.meta["obstypes"])
    log.debug(f"Remove observation types after selection: {' '.join(add_remove_obstypes)}.")
    remove_obstypes.update(add_remove_obstypes)
    _remove_obstype_from_dset(dset, remove_obstypes)
    dset.meta["obstypes"] = selected_obstypes.copy()

    # Remove NaN values of selected observation types
    if config.tech[_SECTION].remove_nan.bool:
        for sys in dset.meta["obstypes"]:
            for type_ in dset.meta["obstypes"][sys]:
                reject_nan = np.isnan(dset.obs[type_][keep_idx])
                if np.any(reject_nan):
                    keep_idx[keep_idx] = np.logical_not(reject_nan)
                    log.debug(f"Remove {np.sum(reject_nan)} NaN values for GNSS '{sys}' observation type {type_}.")
        
    return keep_idx


def _remove_obstype_from_dset(dset: "Dataset", remove_obstypes: Set[str]) -> None:
    """Remove unused observation types from Dataset

    Args:
        dset:             A Dataset containing model data.
        remove_obstypes:  Observations types, which should be removed from Dataset
    """
    for type_ in dset.obs.fields:
        if type_ in remove_obstypes:
            for partition in ["obs", "lli", "snr"]:
                if f"{partition}.{type_}" in dset.fields:
                    del dset[partition][type_]

            for sys in dset.meta["obstypes"]:
                if type_ in dset.meta["obstypes"][sys]:
                    dset.meta["obstypes"][sys].remove(type_)


def _select_observations(obstypes_all, obstypes):
    """Select observations based on GNSS observation priority list

    NOTE: The order how the observation types are saved in 'use_obstypes' is relevant for the further processing, e.g.
          by selection of the code observation type for determination of the satellite transmission time or by
          generation of linear combinations from the observations.

    Args:
        obstypes_all (list):    All observation types defined in RINEX header
        obstypes (dict):        Observation types defined in RINEX header given for each GNSS

    Returns:
        tuple:  with following elements

    =================  ======  ==================================================================================
     Elements           Type    Description
    =================  ======  ==================================================================================
     use_obstypes       dict    Selected observation types for each GNSS related to priority list
     remove_obstypes    set     Set with observation types, which can be removed
    =================  ======  ==================================================================================
    """
    use_obstypes = dict()
    remove_obstypes = set(obstypes_all)
    cfg_freq_type = config.tech.freq_type.str
    cfg_obs_code = config.tech[_SECTION].obs_code.list

    # Convert frequency type in frequency numbers
    try:
        freq_numbers = FREQ_NUMBER_DEF[cfg_freq_type]
    except KeyError:
        log.fatal(f"Configuration option 'freq_type = {cfg_freq_type}' is not valid.")

    # Loop over GNSSs
    for sys in obstypes:
        use_obstypes.update({sys: list()})

        # Loop over observation code
        for obs_code in sorted(cfg_obs_code):
            if obs_code not in OBS_CODE_DEF:
                log.fatal(f"Configuration option 'obs_code= {obs_code}' is not valid.")

            for freq_num in freq_numbers:
                type_ = OBS_CODE_DEF[obs_code] + freq_num
                use_obstypes[sys].append(_select_obstype(sys, type_, obstypes[sys]))

        log.info(f"Selected observation types for GNSS {sys!r}: {', '.join(use_obstypes[sys])}")
        remove_obstypes.difference_update(use_obstypes[sys])

    return use_obstypes, remove_obstypes


def _select_obstype(sys, type_, obstypes):
    """Select GNSS observation type from priority list, which should be used in the Where processing

    Args:
        sys (str):          GNSS.
        type_ (str):        Observation type definition of priority list.
        obstypes (list):    Observation types defined in RINEX header for given GNSS `sys`.

    Returns:
        String:   Selected observation type for given GNSS and priority observation type (L1, C1, L2, C2, L3, C3, ...)
    """

    # TODO: Selection of observations are necessary, for example swisstopo is carrying out following selection:
    #
    # Multi-GNSS analysis at swisstopo with Bernese GNSS Software 5.3
    # Observation code priority list for the assignment on two frequencies
    #
    # Sys Obs   RINEX observation codes and their priority
    # G   L1    L1P L1C L1W L1X
    # G   L2    L2P L2C L2D L2S L2W L2X
    # R   L1    L1P L1C L1X
    # R   L2    L2P L2C L2X
    # E   L1    L1C L1X
    # E   L2    L5Q L5X
    # C   L1    L1I L1X
    # C   L2    L7I L7X
    # G   C1    C1P C1C C1W C1X
    # G   C2    C2P C2C C2D C2S C2W C2X

    # TODO: Priority list has to be improved. After which criteria should the priority be chosen?

    # Loop over observations types defined in priority list
    for obstype in config.where.gnss_obs_priority["{}_{}".format(sys, type_).lower()].list:
        if obstype in obstypes:
            return obstype
