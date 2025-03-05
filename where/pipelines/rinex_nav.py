"""a RINEX navigation files pipeline

Description:
------------
RINEX navigation files are often written separately for each GNSS, whereby the navigation file extension depends on the
GNSS (GPS: *.n, Galileo: *.l, BeiDou: *.c, ...). All these navigation files can be read simultaneously and merged
together. In addition merged (mixed) broadcast ephemeris can be read.

RINEX navigation file manipulation can be used for following tasks:
    - navigation message can be ordered alphabetical after satellite
    - removing of duplicated navigation message records
    - removing of unhealthy satellites
    - removing of defined satellites in configuration file
    - selection of navigation message type e.g. I/NAV or F/NAV Galileo message

TODO: 'report': Information about unhealthy satellites, statistic about number of navigation messages for each
                satellite (also specifically for each navigation message type INAV, FNAV, LNAV)

"""
# Standard library imports
from datetime import datetime
from typing import Dict, List, Tuple, Union

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where import cleaners
from where import writers
from where.lib import config

# The name of this technique
TECH = __name__.split(".")[-1]


@plugins.register_named("options")
def options() -> Tuple[str]:
    """Command line options that can be used to specify this analysis

    Returns:
        Strings specifying command line options.
    """
    return ("--rinex_nav",)


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
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        File variables special for this technique.
    """
    station = config.tech.station.str
    return dict(station=station.lower(), STATION=station.upper())


#
# READ DATA
#
@plugins.register
def read(stage: str, dset: "Dataset") -> None:
    """Read RINEX navigation file

    Args:
        stage:    Name of current stage.
        dset:     A dataset containing the data.
    """
    dset.vars.update(file_vars())
    brdc = apriori.get(
        "orbit",
        rundate=dset.analysis["rundate"],
        system=tuple(config.tech.systems.list),
        station=dset.vars["station"],
        apriori_orbit="broadcast",
        day_offset=config.tech.rinex_nav_day_offset.int,
    )

    # Write either raw or filtered (edit) broadcast ephemeris as Dataset
    if config.tech.get("filter_navigation_message", default=True).bool:
        brdc.dset_raw.write()
        dset.update_from(brdc.dset_edit)
    else:
        dset.update_from(brdc.dset_raw)

    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)


#
# EDIT DATA
#
@plugins.register
def edit(stage: str, dset: "Dataset") -> None:
    """Edit the broadcast ephemeris

    Args:
        stage:      Name of current stage.
        dset:       A dataset containing the data.
    """
    cleaners.apply_removers("removers", dset)

    if util.check_write_level("operational"):
        dset.write_as(stage=stage)


#
# WRITE RESULTS
#
@plugins.register
def write(stage: str, dset: "Dataset") -> None:
    """Write results to file.

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    writers.write(default_dset=dset)

    if util.check_write_level("operational"):
        dset.write_as(stage=stage)
