"""a RINEX observation files pipeline

Description:
------------
RINEX observation file manipulation can be used for following tasks:
    - removing of empty observation type data fields and GNSSs without observations
    - change sampling rate
    - conversion from Android/RINEX2 observation format to RINEX3 observation format
    - writing of GNSS observation report
"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori, cleaners, parsers, postprocessors, writers
from where.lib import config, gnss, log, util

# The name of this pipeline
PIPELINE = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this analysis

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return ("--rinex_obs",)


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

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    sampling_rate = config.tech.sampling_rate.float
    convert_unit = config.tech.convert_unit.bool

    # Read GNSS observation data either from Android raw file or RINEX file
    if config.tech.format.str == "android":
        parser = parsers.parse_key_existing("gnss_android_raw_data", file_vars={**dset.vars, **dset.analysis})
    else:
        version, file_path = gnss.get_rinex_file_version("gnss_rinex_obs")
        log.info(f"Read RINEX file {file_path} with format version {version}.")
        if version.startswith("2"):
            parser = parsers.parse_file(
                "rinex2_obs", file_path=file_path, sampling_rate=sampling_rate, convert_unit=convert_unit
            )
        elif version.startswith("3"):
            parser = parsers.parse_file(
                "rinex3_obs", file_path=file_path, sampling_rate=sampling_rate, convert_unit=convert_unit
            )
        else:
            log.fatal(f"Unknown RINEX format {version} is used in file {file_path}")

    dset.update_from(parser.as_dataset())

    if util.check_write_level("analysis"):                    
        dset.write_as(stage=stage)


#
# ORBIT DETERMINATION
#
@plugins.register
def orbit(stage: str, dset: "Dataset") -> None:
    """Determine GNSS satellite orbits

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    if not config.tech.apriori_orbit.str == "none":
        
        # Remove unused GNSS from observations (Note: This can be necessary for example if no broadcast ephemeris file 
        # is available for a GNSS.)
        cleaners.apply_remover(
                    "gnss_ignore_system", 
                    dset, 
                    systems=set(dset.unique("system")) - set(config.tech.systems.list),
        )
    
        # Calculate satellite position
        orbit = apriori.get(
                "orbit", 
                rundate=dset.analysis["rundate"], 
                system=tuple(dset.unique("system")), 
                station=dset.vars["station"],
        )
        orbit.calculate_orbit(dset, time="time")
    
        # Add satellite position to dataset
        dset.add_posvel("sat_posvel", time=dset.time, system="trs", val=orbit.dset.sat_posvel.trs, other=dset.site_pos)
        
        # Connect site position with satellite orbits needed for determination of elevation and azimuth
        dset.site_pos.other = dset.sat_posvel
        
        #if util.check_write_level("analysis"):
        #    dset.write_as(stage=stage)


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
    # cleaners.apply_editors("editors", dset)
    cleaners.apply_removers("removers", dset)

    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)
    
    
#
# POST-PROCESS DATA
#
@plugins.register
def postprocess(stage, dset):
    """Post-process GNSS data

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    postprocessors.apply_postprocessors("postprocessors", dset)

    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)


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
        dset.write_as(stage="write", dataset_id=0)
