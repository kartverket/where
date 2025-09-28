"""a GNSS SPV pipeline

Description:
------------
Use Doppler measurements from at least four satellites and compute the station velocity and the clock rate.
The algorithm is given by the module spvDoppler class defined in source code rec_velocity_est.py

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.gnss import gnss as mg_gnss

# Where imports
from where.lib import config, gnss, log, util
from where import estimation
from where.models import site, delay
from where import apriori, cleaners, parsers, writers

# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return ("--gnss_vel",)


@plugins.register_named("get_args")
def get_args(rundate, input_args=None):
    """Convert where_runner arguments to where arguments for given date

    Args:
        rundate (date):   The model run date.

    Returns:
        List:   Strings with names of available sessions.
    """
    where_args = set()

    for idx, arg in enumerate(input_args):
        if arg.startswith("--station"):
            stations = arg.split("=")[1].replace(",", " ").split()

            if len(stations) > 1:
                del input_args[idx]
                for station in stations:
                    where_args.add(" ".join(input_args + [f"--station={station}"]))
            else:
                where_args = {" ".join(input_args)}

    return where_args


@plugins.register_named("list_sessions")
def list_sessions(rundate):
    """Get list of sessions for the given rundate

    In the GNSS analysis the sessions correspond to stations. The list of stations is looked up in the config file.

    Args:
        rundate (date):   The model run date.

    Returns:
        List: Strings with names of stations.
    """
    return config.where[TECH].stations.list


@plugins.register_named("file_vars")
def file_vars(file_vars=None):
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """

    # Determine correct interval based on given sampling rate
    sampling_rate = config.tech.sampling_rate.int
    if sampling_rate < 59: # seconds
        interval = str(sampling_rate).zfill(2) + "S"
    elif sampling_rate > 59 & sampling_rate < 3599: # seconds
        interval = str(int(sampling_rate/60)).zfill(2) + "M"
    elif sampling_rate > 3599 & sampling_rate < 86399: # seconds
        interval = str(int(sampling_rate/3600)).zfill(2) + "H"
    elif sampling_rate > 86399 & sampling_rate < 8639999: # seconds
        interval = str(int(sampling_rate/3600/24)).zfill(2) + "D"
    elif sampling_rate > 8639999:
        log.fatal(f"Interval based sampling_rate {sampling_rate} is not defined.")

    return dict(interval=interval, STATION=file_vars["station"].upper())

#
# READ DATA
#
@plugins.register
def read(stage, dset):
    """Read the GNSS RINEX data.

    Following Dataset fields are generated:

   |  Field               | Type              | Description                                                           |
   | :------------------- | :---------------- | :-------------------------------------------------------------------- |
   | <observation type>   | numpy.ndarray     | GNSS observation type data (e.g. C1C, C2W, L1C, L2W, ...) given for   |
   |                      |                   | loss of lock indicator (lli), pseudo-range and carrier phase          |
   |                      |                   | observation (obs) and signal-to-noise-ratio (snr)                     |
   | epoch_flag           | numpy.ndarray     | Epoch flag                                                            |
   | num_satellite_available |  numpy.ndarray | Number of available satellite in each observation epoch               |
   | rcv_clk_offset       | numpy.ndarray     | Receiver clock offset in seconds given for each epoch                 |
   | satellite            | numpy.ndarray     | Satellite PRN number together with GNSS identifier (e.g. G07)         |
   | satnum               | numpy.ndarray     | Satellite PRN number (e.g. 07)                                        |
   | site_pos             | Position          | PositionTable object with given station coordinates either from       |
   |                      |                   | RINEX header or overwritten by 'gnss_station_crd' station coordinate  | 
   |                      |                   | RINEX header)                                                         |
   | station              | numpy.ndarray     | Station name list                                                     |
   | system               | numpy.ndarray     | GNSS identifier                                                       |
   | time                 | Time              | Observation time given as TimeTable object                            |

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    station = dset.vars["station"]
    sampling_rate = config.tech.sampling_rate.float

    # Read GNSS observation data either from Android raw file or RINEX file
    # TODO: Maybe a gnss.py 'obs' modul should be added to ./where/obs?
    if config.tech.format.str == "android":
        parser = parsers.parse_key_existing("gnss_android_raw_data", file_vars={**dset.vars, **dset.analysis})
    else:
        version, file_path = gnss.get_rinex_file_version("gnss_rinex_obs")
        log.info(f"Read RINEX file {file_path} with format version {version}.")
        if version.startswith("2"):
            parser = parsers.parse_file(
                "rinex2_obs",
                file_path=file_path,
                sampling_rate=sampling_rate,
                convert_unit=True,
            )
        elif version.startswith("3"):
            parser = parsers.parse_file(
                "rinex3_obs",
                file_path=file_path,
                sampling_rate=sampling_rate,
                convert_unit=True,
            )
        elif version.startswith("4"):  #TODO: Own RINEX 4 parser should be implemented. But RINEX 3 parser works also.
            parser = parsers.parse_file(
                "rinex3_obs",
                file_path=file_path,
                sampling_rate=sampling_rate,
                convert_unit=True,
            )
        else:
            log.fatal(f"Unknown RINEX format {version} is used in file {file_path}")

    dset.update_from(parser.as_dataset())
    
    # Add site_id and number of satellites per epoch to dataset
    dset.meta.add("site_id", dset.meta["marker_name"].upper(), section=dset.vars["station"])
    dset.add_float(
        "num_satellite_available",
        val=mg_gnss.get_number_of_satellites(dset.system, dset.satellite, dset.time.gps.datetime),
        write_level="operational",
    )

    # Select GNSS observation to process
    cleaners.apply_remover("ignore_epochs", dset)
    cleaners.apply_remover("gnss_select_obs", dset)
    cleaners.apply_remover("gnss_clean_obs", dset)

    # Overwrite station coordinates given in RINEX header
    # TODO: Should be a apriori function with, where a station coordinate can be select for a given station.
    p = parsers.parse_key_existing(parser_name="bernese_crd", file_key="gnss_station_crd")
    sta_crd = p.as_dict()

    if station in sta_crd:
        pos = np.array([sta_crd[station]["pos_x"], sta_crd[station]["pos_y"], sta_crd[station]["pos_z"]])

        # Check station coordinates against RINEX header station coordinates
        limit = 10
        diff = pos - dset.site_pos.trs[0].val
        if not (diff < limit).all():
            log.warn(
                f"Difference between station database (xyz: {pos[0]:.3f} m, {pos[1]:.3f} m, {pos[2]:.3f} m) "
                f"and RINEX header (xyz: {dset.site_pos.trs.x[0]:.3f} m, {dset.site_pos.trs.y[0]:.3f} m, "
                f"{dset.site_pos.trs.z[0]:.3f} m) station coordinates exceeds the limit of {limit} m "
                f"(xyz: {diff[0]:.3f} m, {diff[1]:.3f} m, {diff[2]:.3f} m)."
            )

        # pos = apriori.get("gnss_station_coord", rundate=dset.analysis["rundate"], station=station)
        dset.site_pos[:] = np.repeat(pos[None, :], dset.num_obs, axis=0)
        dset.meta["pos_x"] = sta_crd[station]["pos_x"]
        dset.meta["pos_y"] = sta_crd[station]["pos_y"]
        dset.meta["pos_z"] = sta_crd[station]["pos_z"]

    # Write dataset to file
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage, write_level="analysis")


#
# ORBIT DETERMINATION
#
@plugins.register
def orbit(stage, dset):
    """Determine GNSS satellite orbit

    TODO: Is the workflow for determining the satellite transmission time correct? gLAB determines satellite clock
          correction based on receiver time and not an satellite transmission time. Additionally gLAB does not apply
          relativistic corrections.

    Following Dataset fields are generated:

    | Field                         | Type          | Description                                                    |
    | :---------------------------- | :------------ | :------------------------------------------------------------- |
    | gnss_earth_rotation           | PosVel        | Earth's rotation effect during signal flight time              |
    | navigation_idx                | numpy.ndarray | Indices related to the correct set of broadcast ephemeris for  |
    |                               |               | given observation epochs                                       |
    | sat_clock_bias                | numpy.ndarray | Satellite clock bias in [s]                                    |
    | sat_clock_drift               | numpy.ndarray | Satellite clock drift in [s/s]                                 |
    | sat_clock_drift_rate          | numpy.ndarray | Satellite clock drift rate in [s/s^2]                          |
    | sat_posvel                    | PosVel        | Satellite position and velocity                                |
    | sat_time                      | Time          | Satellite transmission time given as Time object               |
    | used_iode                     | numpy.ndarray | IODE of selected broadcast ephemeris block                     |
    | used_transmission_time        | Time          | Transmission time of selected broadcast ephemeris block        |
    | used_toe                      | Time          | Time of ephemeris (TOE) of selected broadcast ephemeris block  |

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    station = dset.vars["station"]
    orb_flag = config.tech.apriori_orbit.str

    # Second estimate using satellite clock and relativistic clock corrections
    orbit = apriori.get(
        "orbit", rundate=dset.analysis["rundate"], system=tuple(dset.unique("system")), station=station
    )
    if util.check_write_level("analysis"):
        dset_vars = {**dset.vars, **dset.analysis}
        for key in ["label", "stage"]:
            del dset_vars[key]
        orbit.dset_raw.write_as(stage="orbit", label="raw", **dset_vars)
        orbit.dset_edit.write_as(stage="orbit", label="edit", **dset_vars)

    #TODO: Check if it would change anything, if using orbit.calculate_orbit(dset, time="sat_time") instead
    ## First estimate of satellite transmission time
    #sat_time = dset.time - gnss.get_initial_flight_time(dset)

    # Determine initial satellite orbit solution with observation time as approximation
    orbit.calculate_orbit(dset, time="time")

    # Add satellite clock parameter to dataset
    orbit.add_satellite_clock_parameter(dset)

    # TODO: Has it an effect to iterate here to improve satellite transmission time?
    # Determine satellite transmission time based on initial satellite orbit solution
    dset.add_time(
        "sat_time",
        val=dset.time
        - gnss.get_initial_flight_time(
            dset, 
            sat_clock_corr=orbit.dset.gnss_satellite_clock,
            rel_clock_corr=orbit.dset.gnss_relativistic_clock,
        ),
        write_level="operational",
    )

    # Use satellite transmission time for determination of satellite orbits
    orbit.calculate_orbit(dset, time="sat_time")

    # Copy to regular dataset
    dset.add_posvel(
        "sat_posvel", 
        time=dset.sat_time, 
        system="trs", 
        val=orbit.dset.sat_posvel.trs, 
        other=dset.site_pos,
        write_level="operational",
    )
    dset.vars["orbit"] = orb_flag  # Needed e.g. for calling gnss_relativistic_clock model correctly.

    if orb_flag == "broadcast":
        dset.add_float("used_iode", val=orbit.dset.used_iode, write_level="analysis")
        dset.add_time("used_toe", val=orbit.dset.used_toe, write_level="analysis")
        dset.add_time("used_transmission_time", val=orbit.dset.used_transmission_time, write_level="analysis")

    # Connect site position with satellite orbits needed for determination of elevation and azimuth
    dset.site_pos.other = dset.sat_posvel

    # Correct satellite position/velocity due Earth's rotation effect during signal flight time
    sat_pos, sat_vel = gnss.get_earth_rotation(dset)
    log.warn(
        "Correction of satellite position/velocity due to Earth's rotation effect during signal flight time is not"
        " applied."
    )
    # dset.sat_posvel.add_to_itrs(
    #    np.hstack((sat_pos, sat_vel))
    # )  # TODO: Check residuals will be worser by using that. Why?
    dset.add_posvel_delta(
        "gnss_earth_rotation",
        time=dset.sat_time,
        system="trs",
        val=np.hstack((sat_pos, sat_vel)),
        ref_pos=dset.sat_posvel,
        write_level="analysis",
    )
    
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)


#
# EDIT DATA
#
@plugins.register
def edit(stage, dset):
    """Edit GNSS data

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    cleaners.apply_editors("editors", dset)
    cleaners.apply_removers("removers", dset)
    
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)


#
# CALCULATE AND ESTIMATE
#
@plugins.register
def calculate_estimate(stage: str, dset: "Dataset") -> None:
    """Calculate model parameters and estimate

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    max_iterations = config.tech.max_iterations.int

    for iter_num in itertools.count(start=1):

        # MURKS: Workaround: save estimate Dataset before starting a new iteration.
        #       Store estimate results
        if dset.vars["stage"] == "calculate":
            if util.check_write_level("analysis"):
                dset.write_as(stage="estimate", dataset_id=iter_num - 1)

        # CALCULATE
        # -----------
        # Correction of station position in GCRS due to loading and tide effects
        site.calculate_site("site", dset)
        delta_pos = site.add("site", dset)
        dset.site_pos[:] = (dset.site_pos.gcrs + delta_pos[0].gcrs).trs

        # Initialize models given in configuration file by adding model fields to Dataset
        delay.calculate_delay("delay", dset)
        delta_delay = delay.add("delay", dset)

        if "observed" in dset.fields:
            dset.observed[:] = gnss.get_observation(dset)
        else:
            dset.add_float("observed", val=gnss.get_observation(dset), unit="meter/second", write_level="analysis")

        # Get model corrections
        if "calc" in dset.fields:
            dset.calc[:] = delta_delay
        else:
            dset.add_float("calc", val=delta_delay, unit="meter/second", write_level="analysis")

        if "residual" in dset.fields:
            dset.residual[:] = dset.observed - dset.calc
            dset.residual_prefit[:] = dset.residual
        else:
            dset.add_float("residual", val=dset.observed - dset.calc, unit="meter/second", write_level="operational")
            dset.add_float("residual_prefit", val=dset.residual, unit="meter/second", write_level="analysis")

        # Store calculate results
        log.info(f"{dset.num_obs} observations, residual = {dset.rms('residual'):.4f}")
        if util.check_write_level("analysis"):
            dset.write_as(stage="calculate", label=iter_num)

        # ESTIMATE
        # ----------
        # Keep only observation epochs, where 4 satellites are available
        estimation.apply_observation_rejector("gnss_satellite_availability", dset)

        # Get partial derivates
        partial_vectors = estimation.partial_vectors(dset, "estimate_method")

        log.blank()  # Space between iterations for clarity
        log.info(f"Estimating parameters for iteration {iter_num}")
        estimation.call("estimate_method", dset=dset, partial_vectors=partial_vectors, obs_noise=np.ones(dset.num_obs))
        rms = dset.rms("residual")
        log.info(f"{dset.num_obs} observations, postfit residual = {rms:.4f}")

        ## MURKS: Writing of dataset leads to failure in the ongoing processing
        # Store estimate results
        # dset.write_as(stage="estimate", dataset_id=iter_num)
        # dset.read()  # TODO: workaround because caching does not work correctly

        # Detect and remove outliers
        num_obs_before = dset.num_obs
        independent = config.tech.estimate_obs_rejectors_independent.bool
        dset = estimation.apply_observation_rejectors("estimate_obs_rejectors", dset, independent)

        if dset.num_obs == 0:
            log.fatal("No observations available.")
        log.blank()

        if dset.meta["estimate_convergence_status"] and ((num_obs_before - dset.num_obs) == 0):
            log.info(f"Estimation convergence limit of {config.tech.convergence_limit.float:.3e} is fulfilled.")
            break

        if iter_num >= max_iterations:
            break

    # MURKS: better to save dataset after each estimation step
    # Store last estimate results
    if util.check_write_level("analysis"):
        dset.write_as(stage="estimate", label=iter_num)

#
# WRITE RESULTS
#
@plugins.register
def write(stage, dset):
    """Write results to file.

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    writers.write(default_dset=dset)
    if util.check_write_level("operational"):
        dset.write_as(stage="write")
