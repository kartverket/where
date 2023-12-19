"""a Signal-In-Space Range Error (SISRE) pipeline

Description:
------------
SISRE is determined by comparing broadcast against precise orbits and clocks.

"""

# Standard library imports
from datetime import datetime, timedelta
import itertools
from typing import Dict, List, Tuple, Union

# External library imports
import numpy as np

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.math.constant import constant
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where import cleaners
from where import writers
from where.lib import config
from where.lib import gnss
from where.lib import log
from where.lib import util

# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options() -> Tuple[str]:
    """Command line options that can be used to specify this technique

    Returns:
        Strings specifying command line options.
    """
    return ("--sisre",)


@plugins.register_named("get_args")
def get_args(rundate: datetime.date, input_args: Union[None, List[str]]=None) -> List[str]:
    """Convert where_runner arguments to where arguments for given date

    Args:
        rundate:    The model run date.
        input_args: Input arguments.

    Returns:
        Strings with names of available sessions.
    """
    return [" ".join(input_args)] if input_args else list()


@plugins.register_named("file_vars")
def file_vars(file_vars=None) -> Dict[str, str]:
    """File variables that will be available during the running of this technique (via config.files.vars)

    In addition, date and analysis variables are available.

    Returns:
        File variables special for this technique.
    """
    sampling_rate = config.tech.sampling_rate.str
    session_name = config.tech.session_name.str
    station = config.tech.station.str
    return dict(sampling_rate=sampling_rate, session_name=session_name, station=station, STATION=station.upper())


#
# SETUP ORBIT COMPARISON
#
@plugins.register
def setup(stage: str, dset: "Dataset") -> None:
    """Setup Where dataset used for orbit comparison

    The dataset includes the fields 'time', 'system' and 'satellite' like:
        <epoch>  <system>  <satellite>

    Example:

    EPOCH                                  SYSTEM   SATELLITE
    ---------------------------------------------------------------
    datetime.datetime(2016, 3, 1, 1, 0)     'G'      'G01'
    datetime.datetime(2016, 3, 1, 1, 0)     'G'      'G02'
    ...
    datetime.datetime(2016, 3, 1, 1, 0)     'E'      'E01'
    datetime.datetime(2016, 3, 1, 1, 0)     'E'      'E02'
    ...
    datetime.datetime(2016, 3, 1, 1, 5)     'G'      'G01'
    datetime.datetime(2016, 3, 1, 1, 5)     'G'      'G02'
    ...

    The time stamps used for orbit comparison are based on the sampling rate. The 'rundate' at 00:00 o'clock in GPS
    time scale defines the start time. The sampling rate given in the configuration file defines the time offset added
    to start time until 24:00 o'clock in GPS time scale is reached.

    The GNSS satellites are defined in the configuration file, which should be used in the SISRE analysis.

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
   
    
    # Initialize variables
    time = list()
    used_satellites = list()
    used_systems = list()
    dset_satellites = list()
    dset_systems = list()

    dset.vars.update(file_vars())

    # Get options from configuration file
    sampling_rate = config.tech.sampling_rate.float
    systems = config.tech.systems.list
    satellites = config.tech.satellites.list

    # Select satellites for further orbit comparison
    for satellite in satellites:
        if satellite.startswith(tuple(systems)):
            used_satellites.append(satellite)
            used_systems.append(satellite[0])  # Needed for generation of Dataset

    # Prepare generation of Dataset
    hour = 0
    minute = 0
    second = 0
    start_time = datetime(
        dset.analysis["rundate"].year,
        dset.analysis["rundate"].month,
        dset.analysis["rundate"].day,
        hour,
        minute,
        second,
    )
    time_to_save = start_time
    num_satellite = len(used_satellites)
    while time_to_save < start_time + timedelta(days=1):
        time.extend([time_to_save] * num_satellite)
        time_to_save += timedelta(seconds=sampling_rate)
        dset_satellites.extend(used_satellites)  # All satellites has to be given in each epoch
        dset_systems.extend(used_systems)  # All systems identifiers has to be given in each epoch

    # Generate Dataset
    dset.num_obs = len(time)
    dset.add_time("time", val=time, scale="gps", fmt="datetime", write_level="operational")
    dset.add_text("satellite", val=dset_satellites, write_level="operational")
    dset.add_text("system", val=dset_systems, write_level="operational")

    # Get station positions
    # trf = apriori.get('trf', time=dset.time)  # reference_frame='itrf_ssc:2014'
    # stations = trf.sites  ## [k[1] for k, _ in trf.items()]  # TODO: Better solution? trf.stations???
    # if session.upper() in stations:
    #    pass
    # try:
    #   dset.add_position('site_pos', time='time', itrs=trf.pos(('gnss', session.upper())))
    # except KeyError:
    #    pass

    # Write Dataset to file
    if util.check_write_level("analysis"):
        dset.write_as(stage="setup")


#
# EDIT DATA
#
@plugins.register
def edit(stage: str, dset: "Dataset"):
    """Edit the orbit data

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    cleaners.apply_removers("removers", dset)
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)


#
# BROADCAST/PRECISE ORBIT COMPARISON
#
@plugins.register
def calculate(stage: str, dset: "Dataset"):
    """Compare broadcast and precise ephemeris

    Following Dataset fields are generated:

    | Field                  | Type              |  Description
    |------------------------|-------------------|---------------------------------------------------------------------|
    | bias_brdc              | numpy.ndarray     | Satellite bias of broadcast ephemeris in [m]                        |
    | bias_precise           | numpy.ndarray     | Satellite bias of precise orbits in [m]                             |
    | clk_diff               | numpy.ndarray     | Satellite clock correction difference without correction in [m]     |
    | clk_diff_with_dt_mean  | numpy.ndarray     | Satellite clock correction difference corrected for average         |
    |                        |                   | satellite clock offset difference for given GNSS and epoch in [m]   |
    | pco_brdc               | PositionTable     | Phase center offset (PCO) of broadcast ephemeris in [m]             |
    | pco_precise            | PositionTable     | Phase center offset (PCO) of precise orbits in [m]                  |
    | orb_diff               | PositionTable     | Orbit difference in given ITRS [m]                                  |
    | orb_diff_acr           | PositionTable     | Orbit difference given local orbital reference system, that means in|
    |                        |                   | along-track, cross-track and radial (ACR)                           |
    | orb_diff_3d            | numpy.ndarray     | 3D orbit difference based on ACR orbit differences                  |
    | sisre                  | numpy.ndarray     | Signal-in-space range error in [m]                                  |
    | sisre_with_dr_mean     | numpy.ndarray     | Signal-in-space range error with corrected average constellation-   |
    |                        |                   | mean radial orbit error [m]                                         |
    | sisre_orb              | numpy.ndarray     | Orbit-only signal-in-space range error in [m]                       |
    | sisre_orb_with_dr_mean | numpy.ndarray     | Orbit-only signal-in-space range error in with corrected average    |
    |                        |                   | constellation-mean radial orbit error [m]                           |
    | used_iode              | numpy.ndarray     | IODE of selected broadcast ephemeris block                          |
    | used_transmission_time | TimeTable         | Transmission time of selected broadcast ephemeris block             |
    | used_toe               | TimeTable         | Time of ephemeris (TOE) of selected broadcast ephemeris block       |

    The two equations for the determination of the SISRE(orb) and SISRE are given in :cite:`montenbruck2015` (Eq. 1 and
    2 respectively).

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    apply_has_correction = config.tech.apply_has_correction.bool
    sat_clock_related_to = config.tech.sat_clock_related_to.str
    
    ant_precise = apriori.get("gnss_antenna_correction")
    ant_brdc = apriori.get("gnss_antenna_correction", file_key="antex_brdc") #TODO: Do we need something like that? apriori.get("gnss_antenna_correction") if apply_has_correction else apriori.get("gnss_antenna_correction", file_key="antex_brdc")
    brdc, precise = _get_common_brdc_precise_ephemeris(dset)

    # Definition of GNSS frequencies for accessing correct PCOs for SISRE analysis:
    #   for broadcast orbits: PCOs of broadcast orbits are based for Galileo on the average of E1/E5/E6 PCOs defined
    #                         by the European GNSS Service Centre before 29th May 2018. Since 29th May 2018
    #                         (around 4 p.m.) the generation of PCOs for the INAV/FNAV messages was changed. The
    #                         average of the ionosphere-free linear combination of E1/E5a and E1/E5b PCOs is used
    #                         based on PCOs of GSC webpage. The PCOs for GPS are defined by the National Geospatial
    #                         Intelligence Agency (NGA). An ANTEX file is created with these PCOs values, whereby the
    #                         Galileo E1 and respectively GPS L1 is used.
    #   for precise orbits_:  Normally precise Galileo orbits are based on ionosphere-free linear combination E1/E5a
    #                         and for GPS on L1/L2, so that also the PCOs has to be based on that.
    #
    brdc_sys_freq = {"E": "E1", "G": "L1"} #TODO: Do we need something like that? {"E": "E1_E5b", "G": "L1_L2"} if apply_has_correction else {"E": "E1", "G": "L1"} 
    precise_sys_freq = {"E": "E1_E5a", "G": "L1_L2"}

    # Determine satellite position difference (orbit difference in ITRS) related to CoM
    # TODO: If sat_clock_related_to = "apc", then pco_sat.trs should also applied related to APC?
    pco_sat = ant_brdc.satellite_phase_center_offset(brdc.dset, brdc_sys_freq)
    orb_diff = (brdc.dset.sat_posvel.trs - pco_sat.trs) - precise.dset.sat_posvel.trs

    # Determine satellite clock correction difference related to CoM
    bias_brdc, bias_precise = _get_bias(dset, brdc.dset)
    if sat_clock_related_to == "com":
        clk_diff = (brdc.satellite_clock_correction_com(ant_brdc, brdc_sys_freq) - bias_brdc) - (
            precise.satellite_clock_correction_com(ant_precise, precise_sys_freq) - bias_precise
        )
    elif sat_clock_related_to == "apc":
        clk_diff = (brdc.dset.gnss_satellite_clock - bias_brdc) - (precise.dset.gnss_satellite_clock - bias_precise)
    else:
        log.fatal(f"Unknown value {sat_clock_related_to} for option 'sat_clock_related_to'. This option can be " 
                  f"either 'apc' or 'com', which means antenna phase center and respectively center of mass.")
    
    # Apply HAS correction
    if apply_has_correction:
      
        log.info("Apply HAS correction to orbit and clock differences.")
        # TODO: Handling of HAS code bias correction should be improved. Flexible handling of signal type combinations needed.
        dset.add_float(
            "has_code_bias_correction",
            val = _get_bias_has(dset),
            unit="meter",
        )
        orb_diff = orb_diff + dset.has_orbit_correction
        clk_diff = clk_diff + dset.has_clock_correction + dset.has_code_bias_correction

        
    # Calculate SISRE
    dset.add_float("clk_diff", val=clk_diff, unit="meter", write_level="operational")
    dset.add_float("clk_diff_with_dt_mean", val=np.zeros(dset.num_obs), unit="meter", write_level="operational")
    dset.add_float("sisre", val=np.zeros(dset.num_obs), unit="meter", write_level="operational")
    dset.add_float("sisre_with_dr_mean", val=np.zeros(dset.num_obs), unit="meter", write_level="operational")
    dset.add_float("sisre_orb", val=np.zeros(dset.num_obs), unit="meter", write_level="operational")
    dset.add_float("sisre_orb_with_dr_mean", val=np.zeros(dset.num_obs), unit="meter", write_level="operational")
    dset.add_posvel_delta("orb_diff", val=orb_diff, write_level="operational")

    _get_sisre(dset)
    _additional_fields_to_dataset(
        dset, ant_brdc, ant_precise, brdc, precise, bias_brdc, bias_precise, orb_diff, brdc_sys_freq, precise_sys_freq
    )

    # Write raw Dataset to file
    if util.check_write_level("operational"):
        dset.write_as(stage=stage, label="raw")

    # Write cleaned Dataset to file
    outlier_function = config.tech.outlier_function.str
    if outlier_function != "none":
        _reject_outliers(dset, outlier_function)
        if util.check_write_level("operational"):
            dset.write_as(stage=stage, label="clean")


#
# WRITE RESULTS
#
@plugins.register
def write(stage: str, dset: "Dataset"):
    """Write results to file.

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    log.info("Writing model output for SISRE analysis.")
    writers.write(default_dset=dset)


def _additional_fields_to_dataset(
    dset: "Dataset",
    ant_brdc: "AntennaCorrection",
    ant_precise: "AntennaCorrection",
    brdc: "BroadcastOrbit",
    precise: "PreciseOrbit",
    bias_brdc: np.ndarray,
    bias_precise: np.ndarray,
    orb_diff: "PosVelDeltaArray",
    brdc_sys_freq: Dict[str, str],
    precise_sys_freq: Dict[str, str],
):
    """Add additional fields to Dataset.

    Args:
        dset:             Dataset, a dataset containing the data.
        ant_brdc:         Antenna correction object for broadcast orbits
        ant_precise:      Antenna correction object for precise orbits
        brdc:             Broadcast orbit object
        precise:          Precise orbit object
        bias_brdc:        Satellite bias for broadcast orbits
        bias_precise:     Satellite bias for precise orbits
        orb_diff:         Orbit difference
        brdc_sys_freq:    Dictionary with frequency given for GNSS identifier. This is used for selection of correct
                          broadcast orbit ANTEX PCOs.
                              brdc_sys_freq = { <sys_id>: <freq> }
                              (e.g. brdc_sys_freq = {'E': 'E1',  'G': 'L1'} )
        precise_sys_freq: Dictionary with frequency given for GNSS identifier. This is used for selection of correct
                          precise orbit ANTEX PCOs.
                              precise_sys_freq = { <sys_id>: <freq> }
                              (e.g. precise_sys_freq = {'E': 'E1_E5a',  'G': 'L1'} )

    """
    write_level = config.tech.get("write_level", default="operational").str

    # Create Dataset
    fields = {
        "float": {
            "orb_diff_3d": (orb_diff.acr.pos.length, "meter"),
            "used_iode": (brdc.dset.used_iode, ""),
            "diff_trans_toe": (
                ((brdc.dset.used_transmission_time.mjd - brdc.dset.used_toe.mjd) * Unit.day2second),
                "second",
            ),
            "age_of_ephemeris": (((brdc.dset.time.mjd - brdc.dset.used_toe.mjd) * Unit.day2second), "second"),
        },
        "posvel_delta": {
            "pco_precise": (ant_precise.satellite_phase_center_offset(precise.dset, precise_sys_freq), "meter"), # Set precise.dset.meta["pco_sat"]
            "pco_brdc": (ant_brdc.satellite_phase_center_offset(brdc.dset, brdc_sys_freq), "meter"), # Set brdc.dset.meta["pco_sat"]
        },
        "text": {"satellite_type": ant_precise.satellite_type(dset)},
        "time": {"used_transmission_time": brdc.dset.used_transmission_time, "used_toe": brdc.dset.used_toe},
    }
    
    if "used_iode" in dset.fields:
        del fields["float"]["used_iode"]      # Skip adding of "used_iode" to dataset, if it already exists
               
    if write_level == "analysis":
        add_float_fields = {
            "bias_brdc": (bias_brdc, "meter"),
            "bias_precise": (bias_precise, "meter"),
            "clk_brdc_com": (brdc.satellite_clock_correction_com(ant_brdc, brdc_sys_freq), "meter"),
            "clk_precise_com": (precise.satellite_clock_correction_com(ant_precise, precise_sys_freq), "meter"),
        }
        fields["float"].update(add_float_fields)

    for field, value in fields["float"].items():
        if field in dset.fields:
            continue
        dset.add_float(field, val=value[0], unit=value[1], write_level="operational")

    for field, value in fields["posvel_delta"].items():
        if field in dset.fields:
            continue
        
        if field == "pco_precise":
            ref_pos = precise.dset.sat_posvel
        elif field == "pco_brdc":
            ref_pos = brdc.dset.sat_posvel
                        
        dset.add_posvel_delta(
            field, 
            time="time", 
            val=value[0], 
            system="trs", 
            ref_pos=ref_pos, 
            write_level="operational",
        )

    for field, value in fields["text"].items():
        if field in dset.fields:
            continue
        dset.add_text(field, val=value, write_level="operational")

    for field, value in fields["time"].items():
        if field in dset.fields:
            continue
        dset.add_time(field, val=value, write_level="operational")

        #TODO: This does not work: dset._fields[field].write_level = enums.get_value("write_level", "operational")
    

    # Add PCOs for broadcast and precise orbits and seleted configuration settings to Dataset meta
    dset.meta["pco_sat_brdc"] = brdc.dset.meta["pco_sat"]
    dset.meta["pco_sat_precise"] = precise.dset.meta["pco_sat"]
    dset.meta["frequencies"] = config.tech.frequencies.dict
    dset.meta["navigation_message_type"] = config.tech.navigation_message_type.dict
    dset.meta["systems"] = config.tech.systems.list
    dset.meta["service"] = "HAS" if config.tech.apply_has_correction.bool else "OS"
    
    # Remove unnecessary type information
    for type_ in ["frequencies", "navigation_message_type", "systems"]:
        systems = dset.meta[type_] if type(dset.meta[type_]) == list else dset.meta[type_].keys()
        remove_systems = set(systems) - set(dset.unique("system"))
        for sys in remove_systems:
            del dset.meta[type_][sys]


def _get_bias(dset: "Dataset", dset_brdc: "Dataset") -> Tuple[np.ndarray, np.ndarray]:
    """Determine satellite biases for broadcast and precise orbits

    Depending on the GNSS and navigation message type for broadcast and precise ephemeris following DCBs has to be
    applied:

    |System | Type      | Signal   | Broadcast bias | Precise bias                                                      |
    |:------|:----------|:---------|:---------------|:------------------------------------------------------------------|
    |GPS    | LNAV      | G:L1     | tgd            |:math:`-\frac{f^2_{L2}}{f^2_{L1}-f^2_{L2}} DCB^s_{C1W-C2W} + DCB^s_{C1C-C1W}`|
    |       |           | G:L1_L2  | 0              |0                                                                  |
    |Galileo| INAV_E1   | E:E1     | bgd_e1_e5b     |:math:`-\frac{f^2_{E5a}}{f^2_{E1}-f^2_{E5a}} DCB^s_{C1C-C5Q}`      |
    |       | INAV      | E:E1_E5b | 0              |:math:`-\frac{f^2_{E5a}}{f^2_{E1}-f^2_{E5a}} DCB^s_{C1C-C5Q} + \frac{f^2_{E5b}}{f^2_{E1}-f^2_{E5b}} DCB^s_{C1C-C7Q}`|
    |       | FNAV_E5a  | E:E1_E5a | 0              | 0                                                                 |

    Args:
        dset:       Simulated observation data.
        dset_brdc:  Broadcast navigation message data.

    Returns:
        Tuple with following `numpy.ndarray` arrays:

    |Elements       | Description                                                                                      |
    |---------------|--------------------------------------------------------------------------------------------------|
    | bias_brdc     | Satellite bias for broadcast orbits in [m]                                                       |
    | bias_precise  | Satellite bias for precise orbits in [m]                                                         |
    """
    bias_brdc = np.zeros(dset.num_obs)
    bias_precise = np.zeros(dset.num_obs)
    meta_bias_brdc = dict()  # TODO: Exists a better Python solution for this?
    meta_bias_precise = dict()
    for sat in dset.unique("satellite"):
        meta_bias_brdc[sat] = 0
        meta_bias_precise[sat] = 0

    for sys in dset.unique("system"):

        if (sys == "E") and ("E:E1_E5b" in config.tech.frequencies.list):
            dcb = apriori.get("gnss_bias", rundate=dset.analysis["rundate"])
            f_E1 = enums.gnss_freq_E.E1
            f_E5a = enums.gnss_freq_E.E5a
            f_E5b = enums.gnss_freq_E.E5b

            for sat in dset.unique("satellite"):
                if sat.startswith("E"):
                    idx = dset.filter(satellite=sat)
                    dcb_c1c_c5q = (
                        -f_E5a ** 2
                        / (f_E1 ** 2 - f_E5a ** 2)
                        * dcb.get_dsb(sat, "C1C-C5Q", dset.analysis["rundate"])["estimate"]
                    )
                    dcb_c1c_c7q = (
                        f_E5b ** 2
                        / (f_E1 ** 2 - f_E5b ** 2)
                        * dcb.get_dsb(sat, "C1C-C7Q", dset.analysis["rundate"])["estimate"]
                    )
                    bias = (dcb_c1c_c5q + dcb_c1c_c7q) * constant.c
                    bias_precise[idx] = bias
                    meta_bias_precise[sat] = bias

        elif (sys == "E") and ("E:E1" in config.tech.frequencies.list):

            dcb = apriori.get("gnss_bias", rundate=dset.analysis["rundate"])
            f_E1 = enums.gnss_freq_E.E1
            f_E5a = enums.gnss_freq_E.E5a

            for sat in dset.unique("satellite"):
                if sat.startswith("E"):
                    idx = dset.filter(satellite=sat)
                    dcb_c1c_c5q = (
                        -f_E5a ** 2
                        / (f_E1 ** 2 - f_E5a ** 2)
                        * dcb.get_dsb(sat, "C1C-C5Q", dset.analysis["rundate"])["estimate"]
                    )
                    bias = dcb_c1c_c5q * constant.c
                    bias_precise[idx] = bias
                    meta_bias_precise[sat] = bias

                    bias = dset_brdc.bgd_e1_e5b[idx] * constant.c
                    bias_brdc[idx] = bias
                    meta_bias_brdc[sat] = np.mean(bias)

        elif (sys == "G") and ("G:L1" in config.tech.frequencies.list):

            dcb = apriori.get("gnss_bias", rundate=dset.analysis["rundate"])
            f_L1 = enums.gnss_freq_G.L1
            f_L2 = enums.gnss_freq_G.L2

            for sat in dset.unique("satellite"):
                if sat.startswith("G"):
                    idx = dset.filter(satellite=sat)
                    dcb_l1_l2 = (
                        -f_L2 ** 2
                        / (f_L1 ** 2 - f_L2 ** 2)
                        * dcb.get_dsb(sat, "C1W-C2W", dset.analysis["rundate"])["estimate"]
                        + dcb.get_dsb(sat, "C1C-C1W", dset.analysis["rundate"])["estimate"]
                    )
                    bias = dcb_l1_l2 * constant.c
                    bias_precise[idx] = bias
                    meta_bias_precise[sat] = bias

                    bias = dset_brdc.tgd[idx] * constant.c
                    bias_brdc[idx] = bias
                    meta_bias_brdc[sat] = np.mean(bias)

    # Add bias information to Dataset 'meta'
    dset.meta["bias_brdc"] = meta_bias_brdc
    dset.meta["bias_precise"] = meta_bias_precise

    return bias_brdc, bias_precise


def _get_bias_has(dset: "Dataset") -> np.ndarray:
    """Determine Galileo HAS code bias correction needed for correcting Galileo HAS satellite clock 
    correction to be equivalent to precise satellite clocks

    Depending on GNSS frequency type following Galileo HAS code bias corrections has to be applied:

    |System  | Signal   | Galileo HAS code bias corrections |
    |:-------|:---------|:----------------------------------|
    |GPS     | G:L1_L2  | :math:`-\frac{bias_L2 - \frac{f^2_{L1}-f^2_{L2}} bias_L1 + DCB_{C1C,C1P}}{\frac{f^2_{L1}-f^2_{L2}} - 1}` |
    |Galileo | E:E1_E5b | :math:`-\frac{bias_E5b - \frac{f^2_{E1}-f^2_{E5b}} bias_E1}{\frac{f^2_{E1}-f^2_{E5b}} - 1}` |

    Args:
        dset:       Simulated observation data.

    Returns:
        Galileo HAS code bias correction needed for correcting Galileo HAS satellite clock
    """
    correction = np.zeros(dset.num_obs)
    frequencies = config.tech.frequencies.list

    for sys in dset.unique("system"):
        idx = dset.filter(system=sys)

        if (sys == "E") and ("E:E1_E5b" in frequencies):    
            correction[idx] = gnss.ionosphere_free_linear_combination(
                dset.has_code_bias_c1c[idx], 
                dset.has_code_bias_c7q[idx], 
                enums.gnss_freq_E.E1, 
                enums.gnss_freq_E.E5b,
            )
            
        elif (sys == "G") and ("G:L1_L2" in frequencies):
            dcb = apriori.get("gnss_bias", rundate=dset.analysis["rundate"])
            for sat in dset.unique("satellite"):
                if not sys == sat[0:1]:
                    continue
                idx_sat = dset.filter(satellite=sat)
                #TODO: Ideally C1C-CP1 DCB should be used?
                correction[idx_sat] = gnss.ionosphere_free_linear_combination(
                    dset.has_code_bias_c1c[idx_sat] + dcb.get_dsb(sat, "C1C-C1W", dset.analysis["rundate"])["estimate"] * constant.c,
                    dset.has_code_bias_c2p[idx_sat], 
                    enums.gnss_freq_G.L1, 
                    enums.gnss_freq_G.L2,
                )
        else:
            log.fatal(f"Galileo code bias corrections for correcting Galileo HAS clock corrections can not be "
                      f"implemented for frequencies {','.join(frequencies)} (only E:E1_E5b and G:L1_L2 are defined).")

    return correction


def _get_common_brdc_precise_ephemeris(dset: "Dataset") -> Tuple["Dataset", "Dataset"]:
    """Get broadcast and precise orbit data and generate common datasets

    Only common epochs for both broadcast and precise orbit data are used. In addition "observation" dataset is
    adapted.

    Args:
        dset:  A dataset containing the data.        # Apply HAS orbit and clock correction
                #dset.sat_posvel[:] = dset.sat_posvel + dset.has_orbit_correction
                #dset.delay.gnss_satellite_clock[:] = dset.delay.gnss_satellite_clock - dset.has_clock_correction

    Returns:
        Tuple with broadcast and precise apriori data
    """
    
    # Add HAS correction to dataset
    if config.tech.apply_has_correction.bool:
        

        # Add HAS orbit correction. 
        #
        # Note: This step has to be done first to match the broadcast navigation message to the HAS messages and to 
        #       add the "has_gnssiod_orb" field, which is needed to select corresponding broadcast navigation messages
        #       in relation to given reference HAS IOD.
        orbit_has = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            file_key="gnss_has_orb",
            day_offset=0, 
            apriori_orbit="has",
        ) 
        orbit_has.add_orbit_to_dataset(dset)
    
        # Add HAS clock correction
        clock_has = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            file_key="gnss_has_clk",
            day_offset=0, 
            apriori_orbit="has",
        )
        clock_has.add_clock_to_dataset(dset)
        
        # Add HAS code bias correction
        code_bias_has = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            file_key="gnss_has_cb",
            day_offset=0, 
            apriori_orbit="has",
        )
        code_bias_has.add_code_bias_to_dataset(dset)
                  
    # Read, edit and calculate broadcast orbit
    brdc = apriori.get(
        "orbit",
        rundate=dset.analysis["rundate"],
        system=tuple(dset.unique("system")),
        station=dset.vars["station"],
        apriori_orbit="broadcast",
    )
    if util.check_write_level("analysis"):
        brdc.dset_raw.vars = dset.vars.copy()
        brdc.dset_raw.analysis = dset.analysis.copy()
        brdc.dset_raw.write()

        brdc.dset_edit.vars = dset.vars.copy()
        brdc.dset_edit.analysis = dset.analysis.copy()
        brdc.dset_edit.write()
    brdc.calculate_orbit(dset)
    
    # Read, edit and calculate precise orbit
    precise = apriori.get("orbit", rundate=dset.analysis["rundate"], apriori_orbit="precise")
    precise.dset_raw.vars = dset.vars.copy()
    precise.dset_raw.analysis = dset.analysis.copy()
    if util.check_write_level("analysis"):
        precise.dset_raw.vars = dset.vars.copy()
        precise.dset_raw.analysis = dset.analysis.copy()
        precise.dset_raw.write()

        precise.dset_edit.vars = dset.vars.copy()
        precise.dset_edit.analysis = dset.analysis.copy()
        precise.dset_edit.write()
    precise.calculate_orbit(dset)

    # Generate common Dataset for broadcast and precise ephemeris
    # TODO: General solution in dataset.py would be better, then this workaround. Check this solution e.g. on day 12.07.2019.
    keep_idx_dset = np.zeros(dset.num_obs, dtype=bool)
    keep_idx_brdc = np.zeros(brdc.dset.num_obs, dtype=bool)
    keep_idx_precise = np.zeros(precise.dset.num_obs, dtype=bool)

    hash_dset = np.char.add(dset.time.gps.isot, dset.satellite).tolist()
    hash_brdc = np.char.add(brdc.dset.time.gps.isot, brdc.dset.satellite).tolist()
    hash_precise = np.char.add(precise.dset.time.gps.isot, precise.dset.satellite).tolist()
    hash_all = sorted(set(hash_dset) & set(hash_brdc) & set(hash_precise))
    idx_dset = [hash_dset.index(h) for h in hash_all]
    idx_brdc = [hash_brdc.index(h) for h in hash_all]
    idx_precise = [hash_precise.index(h) for h in hash_all]

    keep_idx_dset[idx_dset] = True
    keep_idx_brdc[idx_brdc] = True
    keep_idx_precise[idx_precise] = True
    dset.subset(keep_idx_dset)
    brdc.dset.subset(keep_idx_brdc)
    precise.dset.subset(keep_idx_precise)

    not_common_sat = set(precise.dset.unique("satellite")).difference(set(brdc.dset.unique("satellite")))
    if not_common_sat:
        log.warn(f"The following satellites are not common in brodcast and precise ephemeris: {not_common_sat}")

    # Update dataset vars and analysis variable
    brdc.dset.vars = dset.vars.copy()
    brdc.dset.analysis = dset.analysis.copy()
    precise.dset.vars = dset.vars.copy()
    precise.dset.analysis = dset.analysis.copy()
        
    # Write dataset to file
    if util.check_write_level("analysis"):
        brdc.dset.write_as(stage="calculate", label="brdc")
        precise.dset.write_as(stage="calculate", label="precise")

    return brdc, precise


def _get_sisre(dset: "Dataset") -> None:
    """Calculate SISRE(orb) and SISRE

    Args:
        dset (Dataset):     Model data.
    """
    weight_factor_elev_mask = config.tech.get(key="weight_factor_elev_mask", default=0).int

    # Loop over all GNSS
    for sys in dset.unique("system"):

        # Loop over all observation epochs
        for epoch in dset.unique("time"):
            idx = dset.filter(system=sys, time=epoch)
            da2 = dset.orb_diff.acr.along[idx] ** 2  # Square of along-track difference [m^2]
            dc2 = dset.orb_diff.acr.cross[idx] ** 2  # Square of cross-track difference [m^2]
            dr_mean = np.sum(dset.orb_diff.acr.radial[idx]) / len(
                dset.orb_diff.acr.radial[idx]
            )  # Constellation average
            # radial offset [m]
            dr = dset.orb_diff.acr.radial[idx]
            dt_mean = np.sum(dset.clk_diff[idx]) / len(dset.clk_diff[idx])  # Constellation average clock offset [m]
            dt = dset.clk_diff[idx] - dt_mean  # Satellite clock correction difference [m]

            if weight_factor_elev_mask == 0:
                wr = constant.get("sisre_weight_radial_0deg", source=sys)
                wac2 = constant.get("sisre_weight_along_cross_0deg", source=sys) ** 2
            elif weight_factor_elev_mask == 5:
                wr = constant.get("sisre_weight_radial_5deg", source=sys)
                wac2 = constant.get("sisre_weight_along_cross_5deg", source=sys) ** 2
            else:
                log.fatal(
                    f"Wrong elevation mask with {weight_factor_elev_mask} degree. SISRE weight factors are only defined for 0 or 5 degrees."
                )

            # Orbit-only signal-in-space range error [m]
            dset.sisre_orb[idx] = np.sqrt((wr * dr) ** 2 + wac2 * (da2 + dc2))

            # Signal-in-space range error [m]
            dset.sisre[idx] = np.sqrt((wr * dr - dt) ** 2 + wac2 * (da2 + dc2))

            # Orbit-only signal-in-space range error with constellation average radial offset [m]
            dset.sisre_orb_with_dr_mean[idx] = np.sqrt((wr * (dr - dr_mean)) ** 2 + wac2 * (da2 + dc2))

            # Signal-in-space range error with constellation average radial offset [m]
            dset.sisre_with_dr_mean[idx] = np.sqrt((wr * (dr - dr_mean) - dt) ** 2 + wac2 * (da2 + dc2))

            dset.clk_diff_with_dt_mean[idx] = dt


def _reject_outliers(dset: "Dataset", outlier_function: str) -> None:
    """Detect and remove outliers iteratively with a sigma threshold for each satellite or system based on dataset
       fields

    Outliers are removed from given Dataset 'dset'. The outliers are detected and removed iteratively dependent on
    configuration option 'max_iterations'. An iteration can be done over several Dataset fields for each satellite or
    system. The Dataset fields are defined with option 'outlier_fields'. At the end the SISRE values has to be
    determined again after the outlier rejection. The option 'outlier_function' defines which rejection strategy
    should be used. Following 'outlier_function' exists:
            rms:          Root-mean-square (RMS) for each satellite multiplied by outlier factor used as threshold for
                          removing outliers.
            std:          Mean value for each satellite plus/minus standard deviation multiplied by outlier factor used
                          as threshold for removing outliers.
            std_system:   Mean value for each GNSS plus/minus standard deviation multiplied by outlier factor used as
                          threshold for removing outliers.

    Args:
        dset:     Model data.
    """
    max_iterations = config.tech.get(key="max_iterations", default=3).int + 1
    outlier_factor = config.tech.get(key="outlier_factor", default=3).float
    outlier_fields = config.tech.get(key="outlier_fields", default="clk_diff, sisre_orb").list

    # Detect and remove outliers iteratively
    for iter_num in itertools.count(start=1):

        if iter_num >= max_iterations:
            break

        log.info(f"Detect and remove outliers of SISRE analysis (iteration {iter_num})")

        # Loop over fields
        for field in outlier_fields:

            # Outlier rejection based GNSSs or satellitewise
            if outlier_function == "std_system":
                keep_idx = _reject_outliers_per_system(dset, field, outlier_factor)

            else:
                keep_idx = _reject_outliers_per_satellite(dset, field, outlier_factor, outlier_function)

            dset.subset(keep_idx)
            log.info(f"{sum(~keep_idx)} observations removed")

            # Recalculation of SISRE
            _get_sisre(dset)


def _reject_outliers_per_system(dset: "Dataset", field: str, outlier_factor: float) -> np.ndarray:
    """Detect and remove outliers with a sigma threshold for each GNSS based on dataset fields

    Args:
        dset:            Model data.
        field:           Dataset field.
        outlier_factor:  Factor for definition of sigma-threshold.

    Returns:
        Array with indices. Indices set to False will be rejected.
    """
    # TODO: Can function _reject_outliers_per_satellite and _reject_outliers_per_system be generalized!
    keep_idx = np.ones(dset.num_obs, dtype=bool)

    # Loop over all GNSSs
    for sys in dset.unique("system"):
        sys_idx = dset.filter(system=sys)

        mean = dset.mean(field, idx=sys_idx)
        std = dset.apply(np.std, field, idx=sys_idx)
        field_value = "{:7.4f}  +/-{:7.4f}".format(mean, std)
        keep_sys_idx = (dset[field][sys_idx] > mean - outlier_factor * std) & (
            dset[field][sys_idx] < mean + outlier_factor * std
        )
        outlier_limit = "{:.4f} and smaller than {:.4f}".format(
            mean + outlier_factor * std, mean - outlier_factor * std
        )

        for idx, value in zip(
            np.where(keep_idx >= 0)[0][sys_idx], keep_sys_idx
        ):  # TODO: What would be a better solution for getting indices of an numpy array?
            keep_idx[idx] = value

        log.info(
            f"{sys}: {dset.num_obs - sum(~sys_idx):6d} observations, {field:8s} = {field_value:s}, "
            f"{sum(~keep_sys_idx):6d} rejected observations bigger than {outlier_limit:s}."
        )

        # TODO: Should only be done in logging debug modus. How can it be checked?
        # +DEBUG
        rejected_values = list()
        for time, val in zip(
            dset.time.gps.datetime[sys_idx][np.logical_not(keep_sys_idx)],
            dset[field][sys_idx][np.logical_not(keep_sys_idx)],
        ):
            rejected_values.append(f"{'':40s} {time} {val:>9.4f}")
        log.debug("\n".join(rejected_values))
        # -DEBUG

    return keep_idx


def _reject_outliers_per_satellite(
    dset: "Dataset", field: str, outlier_factor: float, outlier_function: str
) -> np.ndarray:
    """Detect and remove outliers with a sigma threshold for each satellite based on dataset fields

    Args:
        dset:               Model data.
        field:              Dataset field.
        outlier_factor:     Factor for definition of sigma-threshold.
        outlier_function:   Outlier rejection function.

    Returns:
        Array with indices. Indices set to False will be rejected.
    """

    keep_idx = np.ones(dset.num_obs, dtype=bool)

    # Loop over all GNSS satellites
    for sat in dset.unique("satellite"):
        sat_idx = dset.filter(satellite=sat)

        if outlier_function == "rms":
            field_value = "{:7.4f}".format(dset.rms(field, idx=sat_idx))
            keep_sat_idx = np.abs(dset.sisre[sat_idx]) < outlier_factor * dset.rms(field, idx=sat_idx)
            outlier_limit = "{:.4f}".format(outlier_factor * dset.rms(field, idx=sat_idx))
        elif outlier_function == "std":
            mean = dset.mean(field, idx=sat_idx)
            std = dset.apply(np.std, field, idx=sat_idx)
            field_value = "{:7.4f}  +/-{:7.4f}".format(mean, std)
            keep_sat_idx = (dset[field][sat_idx] > mean - outlier_factor * std) & (
                dset[field][sat_idx] < mean + outlier_factor * std
            )
            outlier_limit = "{:.4f} and smaller than {:.4f}".format(
                mean + outlier_factor * std, mean - outlier_factor * std
            )

        for idx, value in zip(
            np.where(keep_idx >= 0)[0][sat_idx], keep_sat_idx
        ):  # TODO: What would be a better solution for getting indices of an numpy array?
            keep_idx[idx] = value

        log.info(
            f"{sat}: {dset.num_obs - sum(~sat_idx):6d} observations, {field:8s} = {field_value:s}, "
            f"{sum(~keep_sat_idx):6d} rejected observations bigger than {outlier_limit:s}."
        )

        # TODO: Should only be done in logging debug modus. How can it be checked?
        # +DEBUG
        rejected_values = list()
        for time, val in zip(
            dset.time.gps.datetime[sat_idx][np.logical_not(keep_sat_idx)],
            dset[field][sat_idx][np.logical_not(keep_sat_idx)],
        ):
            rejected_values.append(f"{'':40s} {time} {val:>9.4f}")
        log.debug("\n".join(rejected_values))
        # -DEBUG

    return keep_idx
