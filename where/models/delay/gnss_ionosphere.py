"""Calculate the delay caused by the ionospheric refraction

Description:
------------
This ionosphere model is applied for single frequency solution, whereby either the Nequick or Klobuchar model can be 
used. Normally the correct ionosphere model is automatically chosen depending on GNSS. But also via configuration can
either the Klobuchar or Nequick model be chosen.

Example:
--------

    from  where.models.delay import gnss_ionosphere
    gnss_ionosphere.gnss_ionosphere(dset)
"""
# External library imports
from enum import Enum
import numpy as np
from typing import Any, Dict, List

# Midgard imports
from midgard.collections.enums import gnss_id_to_3digit_id
from midgard.dev import plugins
from midgard.gnss import gnss, klobuchar
from midgard.parsers import rinex_nav

# Where imports
from where import apriori
from where.lib import config
from where.lib import log

# gl_library imports
from åsgard.glpy import gnss_iono_models


@plugins.register
def gnss_ionosphere(dset: "Dataset") -> np.ndarray:
    """Determine ionosphere correction

    The correct ionospheric model is chosen depending on the GNSS. BeiDou, GPS
    and QZSS uses Klobuchar model and Galileo Nequick model.

    Args:
        dset:   Model data.

    Returns:
        Ionosphere correction in meter
    """
    orbit = apriori.get(
        "orbit", rundate=dset.analysis["rundate"], station=dset.vars["station"], system=tuple(dset.unique("system"))
    )

    corrections = np.zeros(dset.num_obs)
    rundate = dset.analysis["rundate"].strftime("%Y-%m-%d")
    if "iono_para" not in orbit.dset_raw.meta[rundate]:
        log.warn(
            f"Ionospheric parameters do not exists in {' ,'.join(orbit.dset_raw.meta['parser']['file_path'])}. "
            "Ionospheric correction could not be applied."
        )
        return corrections
    iono_para = orbit.dset_raw.meta[rundate]["iono_para"]

    for sys in dset.unique("system"):

        if config.tech.freq_type.str != "single":
            log.warn(f"Ionosphere models can only be used for single-frequency observations.")
            return corrections

        sys_idx = dset.filter(system=sys)
        freq = gnss.obstype_to_freq(sys, dset.meta["obstypes"][sys][0])

        # Use of Nequick model
        if (
            sys == "E" or config.tech.gnss_ionosphere.model.str == "nequick"
        ) and not config.tech.gnss_ionosphere.model.str == "klobuchar":
            log.info("Nequick ionosphere model is used.")
            sys = "E" if sys in ["C", "G", "J"] else sys  # Use Galileo ionosphere parameters for GPS, BeiDou or QZSS
            if config.tech.gnss_ionosphere.filekey_para.str:
                iono_para = _get_iono_para(dset)
                
            nequick_para = iono_para["GAL"]["para"]  # Define input parameter for Nequick model
            corrections[sys_idx] = _nequick(dset, sys_idx, freq, nequick_para)
        elif(
            config.tech.gnss_ionosphere.model.str == "ntcm"
            ):
            log.info("NTCM ionosphere model is used.")
            sys = "E" if sys in ["C", "G", "J"] else sys  # Use Galileo ionosphere parameters for GPS, BeiDou or QZSS
            if config.tech.gnss_ionosphere.filekey_para.str:
                iono_para = _get_iono_para(dset)
                
            ntcm_para = iono_para["GAL"]["para"]  # Define input parameter for Ntcm model
            corrections[sys_idx] = _ntcm(dset, sys_idx, freq, ntcm_para)
        # Use of Klobuchar model
        elif (
            sys in ["C", "G", "J"] or config.tech.gnss_ionosphere.model.str == "klobuchar"
        ) and not config.tech.gnss_ionosphere.model.str == "nequick":
            log.info("Klobuchar ionosphere model is used.")
            sys = "G" if sys == "E" else sys  # Use GPS ionosphere parameters for Galileo
            if config.tech.gnss_ionosphere.filekey_para.str:
                iono_para = _get_iono_para(dset)
                
            iono_alpha = iono_para[gnss_id_to_3digit_id[sys] + "A"]["para"]
            iono_beta = iono_para[gnss_id_to_3digit_id[sys] + "B"]["para"]
            corrections[sys_idx] = _klobuchar(dset, sys, sys_idx, freq, iono_alpha, iono_beta)
            # corrections[sys_idx] = _glpp_klobuchar(dset, sys_idx, freq, iono_alpha, iono_beta)
        elif sys == "I":
            log.fatal("Ionospheric model for IRNSS is not implemented.")
        else:
            log.fatal(f"GNSS '{sys}' is not defined.")

    return corrections


def _get_iono_para(dset: "Dataset") -> Dict[str, Any]:
    """Get ionospheric parameters by reading RINEX navigation file
    
    This is necessary if in used RINEX navigation file needed ionospheric parameters are not given.
    
    Args:
        dset:   Model data.

    Returns:
        Dictionary with ionospheric parameters
    """
    file_path = config.files.path(
            config.tech.gnss_ionosphere.filekey_para.str, file_vars={**dset.vars, **dset.analysis}
    )
    log.debug(f"Read {file_path}")
    parser = rinex_nav.get_rinex2_or_rinex3(file_path)
    return parser.meta["iono_para"]


def _glpp_klobuchar(
                dset: "Dataset", 
                sys_idx: np.ndarray, 
                freq: Enum, 
                iono_alpha: List[float], 
                iono_beta: List[float],
    ) -> np.ndarray:
    """Get Klobuchar ionospheric correction for a Dataset subset based on glpp Klobuchar implementation
    
    Args:
        dset:           Model data.
        sys_idx:        Index mask array for selecting Dataset observation for given GNSS.
        freq:           Frequency in Hz.
        iono_alpha:     Klobuchar ionospheric coefficients.
        iono_beta:      Klobuchar ionospheric coefficients.

    Returns:
        Ionospheric corrections
    """

    # Initialze gnss_iono_models object
    klob = gnss_iono_models.Klobuchar(error_devices=("console",))

    # Define input parameter for Klobuchar model
    rec_pos = dset.site_pos.itrs.mean(axis=0)
    params = {
        "alpha0": iono_alpha[0],
        "alpha1": iono_alpha[1],
        "alpha2": iono_alpha[2],
        "alpha3": iono_alpha[3],
        "beta0": iono_beta[0],
        "beta1": iono_beta[1],
        "beta2": iono_beta[2],
        "beta3": iono_beta[3],
        "a_sub0": 0,
        "a_sub1": 0,
        "a_sub_ot": 0,
        "delta_tls": 0,
        "delta_tlsf": 0,
        "ion_time": 0,
        "wn_sub_t": 0,
        "wn_sub_lsf": 0,
        "dn": 0,
    }

    # .. Get ionospheric delay of Klobuchar model
    corr_subarray = np.zeros(len(dset.time.gps_ws.seconds[sys_idx]))
    for idx, (sat_pos, time) in enumerate(zip(dset.sat_posvel.itrs_pos[sys_idx], dset.time.gps.gps_seconds[sys_idx])):
        corr_subarray[idx] = klob.get_iono_delay(rec_pos, sat_pos, time, freq, params)

    return corr_subarray


def _nequick(
            dset: "Dataset", 
            sys_idx: np.ndarray, 
            freq: Enum, 
            nequick_para,
    ) -> np.ndarray:
    """Get Nequick ionospheric correction for a Dataset subset
    
    Args:
        dset:           Model data.
        sys_idx:        Index mask array for selecting Dataset observation for given GNSS.
        freq:           Frequency in Hz.
        nequick_para:   Nequick ionospheric parameters.

    Returns:
        Ionospheric corrections
    """
    
    # Initialze gnss_iono_models object
    nq = gnss_iono_models.NeQuickG(error_devices=("console",))
    
    rec_pos = dset.site_pos.trs.val.mean(axis=0)  # MURKS: Is that correct?
    
    params = {
        "week_number": int(dset.time.gps_ws.week.mean()),
        "time_of_week": int(dset.time.gps_ws.seconds.mean()),
        "ai0": nequick_para[0],
        "ai1": nequick_para[1],
        "ai2": nequick_para[2],
        "sf": int(nequick_para[3]),
    }

    # .. Get ionospheric delay of Nequick model
    corr_subarray = np.zeros(len(dset.time.gps_ws.seconds[sys_idx]))
    for idx, (prn, sat_pos, time) in enumerate(
        zip(dset.satnum[sys_idx], dset.sat_posvel.trs.pos.val[sys_idx], dset.time.gps.gps_seconds[sys_idx])
    ):
        params.update({"prn": int(prn)})
        corr_subarray[idx] = nq.get_iono_delay(rec_pos, sat_pos, time, freq, params)
        
    return corr_subarray

def _ntcm(
            dset: "Dataset", 
            sys_idx: np.ndarray, 
            freq: Enum, 
            ntcm_para,
    ) -> np.ndarray:
    """Get NTCM ionospheric correction for a Dataset subset
    
    Args:
        dset:           Model data.
        sys_idx:        Index mask array for selecting Dataset observation for given GNSS.
        freq:           Frequency in Hz.
        ntcm_para:      NTCM ionospheric parameters.

    Returns:
        Ionospheric corrections
    """
    
    # Initialze gnss_iono_models object
    ntcm = gnss_iono_models.NTCMG(error_devices=("console",))
    
    rec_pos = dset.site_pos.trs.val.mean(axis=0)  # MURKS: Is that correct?
    
    params = {
        "week_number": int(dset.time.gps_ws.week.mean()),
        "time_of_week": int(dset.time.gps_ws.seconds.mean()),
        "ai0": ntcm_para[0],
        "ai1": ntcm_para[1],
        "ai2": ntcm_para[2],
        "sf": int(ntcm_para[3]),
    }

    # .. Get ionospheric delay of NTCM model
    corr_subarray = np.zeros(len(dset.time.gps_ws.seconds[sys_idx]))
    for idx, (prn, sat_pos, time) in enumerate(
        zip(dset.satnum[sys_idx], dset.sat_posvel.trs.pos.val[sys_idx], dset.time.gps.gps_seconds[sys_idx])
    ):
        params.update({"prn": int(prn)})
        corr_subarray[idx] = ntcm.get_iono_delay(rec_pos, sat_pos, time, freq, params)
        
    return corr_subarray

def _klobuchar(
            dset: "Dataset", 
            sys: str,
            sys_idx: np.ndarray, 
            freq: Enum, 
            iono_alpha: List[float], 
            iono_beta: List[float],
    ) -> np.ndarray:
    
    
    """Get Klobuchar ionospheric correction for a Dataset subset
    
    Args:
        dset:           Model data.
        sys:            GNSS identifier.
        sys_idx:        Index mask array for selecting Dataset observation for given GNSS.
        freq:           Frequency in Hz.
        iono_alpha:     Klobuchar ionospheric coefficients.
        iono_beta:      Klobuchar ionospheric coefficients.

    Returns:
        Ionospheric corrections
    """

    # Define input parameter for Klobuchar model
    rec_llh_pos = dset.site_pos.llh.mean(axis=0)
    freq_l1 = gnss.obstype_to_freq(sys, "L1")  # TODO: This should be done by klobuchar routine.

    # Get ionospheric delay of Klobuchar model
    corr_subarray = np.zeros(len(dset.time.gps_ws.seconds[sys_idx]))
    for idx, (gpssec, az, el) in enumerate(
        zip(dset.time.gps_ws.seconds[sys_idx], dset.site_pos.azimuth[sys_idx], dset.site_pos.elevation[sys_idx])
    ):
        try:
            delay, _ = klobuchar.klobuchar(gpssec, iono_alpha + iono_beta, rec_llh_pos, az, el, freq_l1, freq)
        except ValueError as err:
            log.warn(f"Problem with ionosphere model: {err}")
        corr_subarray[idx] = delay

    return corr_subarray
