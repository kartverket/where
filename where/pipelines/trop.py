"""a Troposphere pipeline

Description:
------------

TODO

"""
import itertools

# External library imports
import numpy as np
from midgard.math import interpolation

# Where imports
from where import apriori
from where import obs
from where import data
from where import cleaners
from where import estimation
from where import parsers
from where.lib import config
from where.lib import exceptions
from where.lib import files
from where.lib import log
from where.lib import plugins
from where.lib import util
from where import models
from where import writers

# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options():
    """Command line options that can be used to specify this technique

    Returns:
        Tuple:  Strings specifying command line options.
    """
    return "-t", "--trop"


@plugins.register_named("list_sessions")
def list_sessions(rundate):
    """Sessions available for the given rundate

    Args:
        rundate (date):   The model run date.

    Returns:
        List:   Strings with names of available sessions.
    """
    return config.where[TECH].stations.list


@plugins.register_named("validate_session")
def validate_session(rundate, session):

    return session


@plugins.register_named("file_vars")
def file_vars():
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        Dict:  File variables special for this technique.
    """
    return dict()


@plugins.register
def read(rundate, session, prev_stage, stage):
    """Read VLBI data

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.

    Returns:
        Bool: True if data are available for the session, False otherwise
    """
    # Parse SINEX TRO file
    trop_data = parsers.parse_key("trop_files")
    sessiondata = trop_data.as_dict()[session]
    # Store coordinates
    pos = sessiondata.pop("pos")
    frame = sessiondata.pop("trf_system")
    dset = data.Dataset.anonymous(num_obs=len(sessiondata["epoch"]))
    dset.add_time("time", val=sessiondata.pop("epoch"), scale="gps", format="datetime")
    for field, value in sessiondata.items():
        dset.add_float(field, val=value)

    dset.write_as(rundate=rundate, tech=TECH, stage=stage, dataset_name=session, dataset_id=0)
    log.info(f"Parsed {dset.num_obs} observations")
    log.info(f"Reading ZTD data for station {session}")
    return True


@plugins.register
def calculate(rundate, session, prev_stage, stage):
    """Estimate model parameters

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """
    # import IPython
    # IPython.embed()
    sampling_rate = config.tech.sampling_rate.int
    interp_type = config.tech.interpolator.str
    dset = data.Dataset(rundate, tech=TECH, stage=prev_stage, dataset_name=session, dataset_id="last")
    sec = dset.time.sec_of_day
    newtime = np.arange(min(sec), max(sec) + 1, sampling_rate)
    newdataset = data.Dataset.anonymous(len(newtime))
    day0 = dset.time.mjd_int[0]
    newdataset.add_time("time", val=day0, val2=newtime / 86400, scale=dset.time.scale, format="mjd")
    for field in dset.fields:
        if field == "time":
            continue
        newval = interpolation.interpolate(sec, dset[field], newtime, kind=interp_type)
        newdataset.add_float(field, val=newval)
    log.info(f"Interpolating data for station {session} at {sampling_rate} s rate")
    # Store results
    newdataset.write_as(rundate=rundate, tech=TECH, stage=stage, dataset_name=session, dataset_id=0)


@plugins.register
def write_result(rundate, session, prev_stage, stage):
    """Write results to file

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        rundate (Datetime):  The model run date.
        session (String):    Name of session.
        prev_stage (String): Name of previous stage.
        stage (String):      Name of current stage.
    """
    writers.write(default_stage=prev_stage)
