"""a RINEX observation files pipeline

Description:
------------
RINEX observation file manipulation can be used for following tasks:
    - removing of empty observation type data fields and GNSSs without observations
    - change sampling rate
    - conversion from Android/RINEX2 observation format to RINEX3 observation format

TODO: 'report': GNSS quality information

"""

# Standard library imports
from datetime import datetime

# WHERE imports
from where import apriori
from where import cleaners
from where import data
from where import parsers
from where import writers
from where.lib import config
from where.lib import gnss
from where.lib import log
from where.lib import plugins
from where.reports import report

# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this analysis

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return ("--rinex_obs",)


@plugins.register_named("file_vars")
def file_vars():
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    station = config.tech.station.str
    return dict(station=station, STATION=station.upper())


#
# READ DATA
#
@plugins.register
def read_data(rundate, session, prev_stage, stage):
    """Read the GNSS RINEX data.

    Args:
        rundate (datetime.datetime): The model run date.
        session (str):               Name of session.
        prev_stage (str):            Name of previous stage.
        stage (str):                 Name of current stage.
    """
    station = config.tech.station.str
    log.info("Parsing {}-data from {}", station, rundate.strftime(config.FMT_date))

    # Read RINEX file in format version 2 or 3 or raw data from Android
    file_vars = dict(station=station, STATION=station.upper())

    # TODO: Better solution for handling android data
    if config.tech.format.str == "android":
        parser = parsers.parse("gnss_android_raw_data", rundate=rundate, station=station)
    else:
        version, filepath = gnss.get_rinex_file_version("gnss_rinex_obs", file_vars)
        log.info("RINEX file format {}.", version)
        if version.startswith("2"):
            parser = parsers.parse("rinex2_obs", rundate=rundate, station=station)
        elif version.startswith("3"):
            parser = parsers.parse("rinex3_obs", rundate=rundate, station=station)
        else:
            log.fatal("Unknown RINEX format {} is used in file {}.", version, filepath)

    dset = data.Dataset(
        rundate, **dict(parser.vars, tech=TECH, stage=stage, dataset_name="", dataset_id=0, empty=True)
    )
    parser.write_to_dataset(dset)
    dset.write()


#
# EDIT DATA
#
@plugins.register
def edit_data(rundate, session, prev_stage, stage):
    """Edit GNSS data

    Args:
        rundate (datetime.datetime): The model run date.
        session (str):               Name of session.
        prev_stage (str):            Name of previous stage.
        stage (str):                 Name of current stage.
    """
    station = config.tech.station.str
    dset = data.Dataset(rundate, tech=TECH, stage=prev_stage, dataset_name="", dataset_id="last")
    cleaners.apply_removers("removers", dset)
    if dset.num_obs == 0:
        log.warn("No observations are available for station '{}'.", station)
        return False

    dset.write_as(stage=stage)


#
# WRITE RESULTS
#
@plugins.register
def write_result(rundate, session, prev_stage, stage):
    """Write results to file.

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        rundate (datetime.datetime): The model run date.
        session (str):               Name of session.
        prev_stage (str):            Name of previous stage.
        stage (str):                 Name of current stage.
    """
    log.info("Writing model output for station '{}'", config.tech.station.str)
    writers.write(default_stage=prev_stage)
