"""a GNSS High Accuracy Service (Galileo) pipeline

Description:
------------

TODO

"""

# Standard library imports
from datetime import datetime

# External library imports
from typing import Dict, List, Tuple, Union

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where import cleaners
from where import writers
from where.lib import config
from where.lib import util

# The name of this technique
PIPELINE = __name__.split(".")[-1]

TEST = False


@plugins.register_named("options")
def options() -> Tuple[str]:
    """Command line options that can be used to specify this technique

    Returns:
        Strings specifying command line options.
    """
    return "--gnss_has",


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
    """Read the Galileo HAS data.

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    dset.vars.update(file_vars())  
    has = apriori.get(
        "orbit", 
        rundate=dset.analysis["rundate"], 
        file_key=dset.vars["file_key"],
        day_offset=config.tech.has_day_offset.int, 
        apriori_orbit="has",
    )
    
    # Write either raw or cleaned (edit) HAS messages as Dataset
    if config.tech.clean_has_message.bool:
        has.dset_raw.write()
        dset.update_from(has.dset_edit)
    else:
        dset.update_from(has.dset_raw)
    
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)
    

#
# EDIT DATA
#
@plugins.register
def edit(stage: str, dset: "Dataset") -> None:
    """Edit the Galileo HAS data.

    Args:
        stage:      Name of current stage.
        dset:       A dataset containing the data.
    """
    cleaners.apply_removers("removers", dset)
    if util.check_write_level("analysis"):
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
        dset.write_as(stage="write")
    
    

