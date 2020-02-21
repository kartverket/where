"""a GNSS SPV pipeline

Description:
------------
Use Doppler measurements from at least four satellites and compute the station velocity and the clock rate.
The algorithm is given by the module spvDoppler class defined in source code rec_velocity_est.py

"""

# Standard library imports
import time

# External library imports
import numpy as np
from math import pi

# Midgard imports
from midgard.dev import plugins
from midgard.gnss.compute_dops import compute_dops

# Where imports
from where.lib import config
from where.lib import gnss
from where.lib import log
from where.models import site, delay
from where import apriori
from where import cleaners
from where import parsers
from where import writers


# ========================================================
# Midgard imports SPV class
# ========================================================
from midgard.gnss.rec_velocity_est import spvDoppler

# coordinate transformation
# import navpy


# The name of this technique
TECH = __name__.split(".")[-1]

TEST = False


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return ("--gnss_spv",)


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
            del input_args[idx]
            stations = arg.split("=")[1].replace(",", " ").split()

            if len(stations) > 1:
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
def file_vars():
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    if not config.analysis.station.str:
        log.fatal(
            f"Station name is not defined. Use 'station' option for defining which station should be used in {TECH.upper()} analysis."
        )
    station = config.analysis.station.str
    return dict(station=station, STATION=station.upper())


# ====================================================== #
# STEP1: read data based on configuration  parameters	 #
# ====================================================== #
@plugins.register
def read(stage, dset):
    """Read the GNSS RINEX data.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    dset.vars.update(file_vars())
    station = dset.vars["station"]
    sampling_rate = config.tech.sampling_rate.float

    # Read GNSS observation data either from Android raw file or RINEX file
    # TODO: Maybe a gnss.py 'obs' modul should be added to ./where/obs?
    if config.tech.format.str == "android":
        parser = parsers.parse_key("gnss_android_raw_data", rundate=dset.analysis["rundate"], station=station)
    else:
        version, file_path = gnss.get_rinex_file_version("gnss_rinex_obs")
        log.info(f"Read RINEX file {file_path} with format version {version}.")
        if version.startswith("2"):
            parser = parsers.parse_key("rinex2_obs", file_path=file_path, sampling_rate=sampling_rate)
        elif version.startswith("3"):
            parser = parsers.parse_file("rinex3_obs", file_path=file_path, sampling_rate=sampling_rate)
        else:
            log.fatal(f"Unknown RINEX format {version} is used in file {file_path}")

    dset.update_from(parser.as_dataset())

    # Select GNSS observation to process
    cleaners.apply_remover("gnss_select_obs", dset)

    # Overwrite station coordinates given in RINEX header
    # TODO: Should be a apriori function with, where a station coordinate can be select for a given station.
    #      "check_coordinate"/"limit" -> station coordinate given in RINEX header and "database" could be checked
    #                                 -> warning could be given

    p = parsers.parse_key(parser_name="gnss_bernese_crd", file_key="gnss_station_crd")
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
    dset.write_as(stage=stage)


# ============================================================ #
# STEP2: Orbit determination (satellite position and velocity  #
# ============================================================ #
#
# ORBIT DETERMINATION
#
@plugins.register
def orbit(stage, dset):
    """Determine GNSS satellite orbit

    TODO: Is the workflow for determining the satellite transmission time correct? gLAB determines satellite clock
          correction based on receiver time and not an satellite transmission time. Additionally gLAB does not apply
          relativistic corrections.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    station = dset.vars["station"]
    orb_flag = config.tech.apriori_orbit.str

    # First estimate of satellite transmission time
    sat_time = dset.time - gnss.get_initial_flight_time(dset)

    # Second estimate using satellite clock and relativistic clock corrections
    orbit = apriori.get(
        "orbit", rundate=dset.analysis["rundate"], system=tuple(dset.unique("system")), station=station
    )
    orbit.dset_raw.write_as(stage=stage, station=station, label="raw")
    orbit.dset_edit.write_as(stage=stage, station=station, label="edit")

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
            dset, sat_clock_corr=orbit.dset.gnss_satellite_clock, rel_clock_corr=orbit.dset.gnss_relativistic_clock
        ),
    )

    # Use satellite transmission time for determination of satellite orbits
    orbit.calculate_orbit(dset, time="sat_time")

    # Copy to regular dataset
    dset.add_posvel("sat_posvel", time=dset.sat_time, system="trs", val=orbit.dset.sat_posvel.trs, other=dset.site_pos)
    dset.add_float("delay.gnss_satellite_clock", val=-orbit.dset.gnss_satellite_clock, unit="meter")
    dset.add_float("delay.gnss_relativistic_clock", val=-orbit.dset.gnss_relativistic_clock, unit="meter")
    dset.vars["orbit"] = orb_flag  # Needed e.g. for calling gnss_relativistic_clock model correctly.

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
    )
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
    dset.write_as(stage=stage)


#
# CALCULATE
#
@plugins.register
def calculate(stage, dset):

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
        dset.observed[:] = gnss.get_code_observation(dset)
    else:
        dset.add_float("observed", val=gnss.get_code_observation(dset), unit="meter")

    # Get model corrections
    if "calc" in dset.fields:
        dset.calc[:] = delta_delay
    else:
        dset.add_float("calc", val=delta_delay, unit="meter", write_level="operational")

    if "residual" in dset.fields:
        dset.residual[:] = dset.observed - dset.calc
    else:
        dset.add_float("residual", val=dset.observed - dset.calc, unit="meter")

    # Store calculate results
    log.info(f"{dset.num_obs} observations, residual = {dset.rms('residual'):.4f}")
    dset.write_as(stage="calculate", dataset_id=0)


# ===========================================
# 	Single Point Velocity by Doppler
#   		(ECEF = ITRS)
# ===========================================
@plugins.register
def spv_doppler(stage, dset):
    """Calculate model parameters and estimate

    Args:
        stage (str)   :       Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """

    logger = print
    dop_obs = np.zeros(50)
    n = dop_obs.size
    H = np.zeros(shape=(n, 4))
    x = dx = np.zeros(4)
    Qx = np.zeros(shape=(4, 4))

    # ================================================================
    # 1. create an object of type spvDoppler and display the content
    #    of the attributes - Constructor (Creator and Initializer)
    # ================================================================
    spv_obj = spvDoppler(dop_obs, H, x, dx, Qx)
    spv_obj.display_attrb()

    # ================================================================
    # 2. retrieve input data from Where dataset dset object
    # ================================================================
    rcv_pos = np.array([dset.meta["pos_x"], dset.meta["pos_y"], dset.meta["pos_z"]])
    n_items = dset.num_obs
    spv_sol = np.zeros(shape=(n_items, 4))
    spv_vel = np.zeros((n_items))

    # administrate the output results from SPV
    n_out_items_cnt = 20
    spv_tmp_vel = np.zeros(n_out_items_cnt)
    gnss_spv_results = np.zeros(shape=(n_items, n_out_items_cnt))
    MJD = sorted(set(dset.time.gps.mjd))  # get the observed time
    for epoch in MJD:
        idx = dset.time.gps.mjd == epoch
        valid_sats = dset.satellite[idx]
        sat_cnt = valid_sats.size
        Az = dset.site_pos.azimuth[idx]  # Az in [0, 2*pi]
        Az[Az < .0] += 2 * pi
        El = dset.site_pos.elevation[idx]
        sat_vel = dset.sat_posvel.trs.vel[idx,]
        sat_dts = dset.sat_clock_drift[idx]
        sat_pos = dset.sat_posvel.trs.pos[idx,]  # ECEF
        dop_obs = gnss.get_code_observation(dset)[idx]

        # ================================================================
        # 3. Estimate the receiver velocity based on Doppler measurements
        # ================================================================
        if sat_cnt > 3:

            spv_obj.est_spv_by_doppler(dop_obs, sat_pos, sat_vel, sat_dts, rcv_pos, Az, El, valid_sats, logger)

            # ========================================================
            #  A.  Update the solution and the dataset dset
            # ========================================================
            xyz_vel = spv_obj.x[0:3]
            est_rcv_vel = np.linalg.norm(xyz_vel)
            h_est_rcv_vel = np.linalg.norm(spv_obj.x[0:2])

            spv_vel[idx] = est_rcv_vel
            spv_sol[idx,] = spv_obj.x

            # ===============================
            # Compute the geometry factors
            # ===============================
            gdop, pdop, tdop, hdop, vdop = compute_dops(Az, El)

            # get the transformation matrix [XYZ ---> ENU]
            mat_T = spv_obj.rot_mat_T
            enu_vel = xyz_vel @ mat_T

            b_flag = False
            if b_flag:
                logger(
                    f" gnss_SPV():: Geometry factors:: GDOP={gdop}, PDOP={pdop}, HDOP={hdop}, VDOP={vdop}, --> [Epoch={epoch},  observed={sat_cnt}]"
                )
                logger(f"\n\n gnss_SPV():: mat_T:\n")
                print(mat_T)
                logger(f"\n\n gnss_SPV():: enu_vel:\t")
                print(ENU_vel)

            logger(
                f"\t\t\t   gnss_SPV():: computed receiver velocity by Doppler measurements(3D)::{est_rcv_vel} [m/s], 2D::{h_est_rcv_vel} [m/s]"
            )
            time.sleep(1)

            # ==============================================================================
            # Update dataset with results generated by spv_doppler process.
            # This includes five groups and are defined as follow:
            #   GR-I :  DATE, MJD, GPS-WEEK, and GPSSEC
            #   GR-II:  SOL_VAL, satCNT, 3D and 2D velocities
            #   GR-III: receiver velocity in XYZ and ENU coordinates
            #   GR-IV:  Geometry factors(GDOP, PDOP, HDOP, and VDOP)
            #   GR-V:   Solution VCM (COV_XX, COV_XY, COV_XZ, COV_YY, COV_YZ and COV_ZZ)
            # ==============================================================================
            # Handle GR-V: Solusion variance-covariance matrix elements
            prec_mat = spv_obj.Qx
            cov_xx = prec_mat[0, 0]
            cov_yy = prec_mat[1, 1]
            cov_zz = prec_mat[2, 2]
            # variances
            cov_xy = prec_mat[0, 1]
            cov_xz = prec_mat[0, 2]
            cov_yz = prec_mat[1, 3]
            # covariances

            # validation test passes if Chi-square and Geometry factos are tested
            tmp_var = 1
            if spv_obj.sol_val:
                tmp_var = 1
            else:
                tmp_var = 0

            if pdop > config.tech.gnss_pdop.pdop_limit.float:
                tmp_var = 0

            # Update SPV data
            spv_tmp_vel = [
                tmp_var,
                sat_cnt,
                est_rcv_vel,
                h_est_rcv_vel,
                xyz_vel[0],
                xyz_vel[1],
                xyz_vel[2],
                enu_vel[0],
                enu_vel[1],
                enu_vel[2],
                gdop,
                pdop,
                hdop,
                vdop,
                cov_xx,
                cov_yy,
                cov_zz,
                cov_xy,
                cov_xz,
                cov_yz,
            ]
            gnss_spv_results[idx,] = spv_tmp_vel

        else:
            logger(f" gnss_SPV():: Needed at least 4 observations, observed={sat_cnt} --> Epoch={epoch}")

    # =====================================================
    # 4. We are done, update dataset with gnss_spv_results
    # =====================================================

    # used for test
    dset.add_float("gnss_spv_vel", val=spv_vel)

    # Handle  GR-II= 4 items
    dset.add_float("sol_val", val=gnss_spv_results[:, 0])
    dset.add_float("num_satellite_used", val=gnss_spv_results[:, 1])
    dset.add_float("3d_vel", val=gnss_spv_results[:, 2], unit="meter/second")
    dset.add_float("2d_vel", val=gnss_spv_results[:, 3], unit="meter/second")

    # Handle  GR-III= 6 items
    dset.add_float("x_vel", val=gnss_spv_results[:, 4], unit="meter/second")
    dset.add_float("y_vel", val=gnss_spv_results[:, 5], unit="meter/second")
    dset.add_float("z_vel", val=gnss_spv_results[:, 6], unit="meter/second")
    dset.add_float("e_vel", val=gnss_spv_results[:, 7], unit="meter/second")
    dset.add_float("n_vel", val=gnss_spv_results[:, 8], unit="meter/second")
    dset.add_float("u_vel", val=gnss_spv_results[:, 9], unit="meter/second")

    # Handle  GR-IV = 4 items
    dset.add_float("gdop", val=gnss_spv_results[:, 10], unit="meter")
    dset.add_float("pdop", val=gnss_spv_results[:, 11], unit="meter")
    dset.add_float("hdop", val=gnss_spv_results[:, 12], unit="meter")
    dset.add_float("vdop", val=gnss_spv_results[:, 13], unit="meter")

    #  Handle  GR-V: 6 items
    dset.add_float("c_xx", val=gnss_spv_results[:, 14], unit="meter")
    dset.add_float("c_yy", val=gnss_spv_results[:, 15], unit="meter")
    dset.add_float("c_zz", val=gnss_spv_results[:, 16], unit="meter")
    dset.add_float("c_xy", val=gnss_spv_results[:, 17], unit="meter")
    dset.add_float("c_xz", val=gnss_spv_results[:, 18], unit="meter")
    dset.add_float("c_yz", val=gnss_spv_results[:, 19], unit="meter")

    dset.write_as(stage="spv_doppler")

    # compute some statistics for test !
    my_min = np.min(spv_vel)
    my_max = np.max(spv_vel)

    log.debug(f" gnss_SPV():: Statistics MIN= {my_min}  MAX= {my_max}")


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
