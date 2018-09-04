"""Get current master file for VLBI

Description:
------------


"""
# Where imports
from where import parsers
from where.lib import config
from where.lib import log
from where.lib import plugins


@plugins.register
def get_master_file(rundate=None):
    """Read master file

    If rundate is not specified, used file_vars that are already set.
    """
    file_vars = None if rundate is None else config.date_vars(rundate)
    return parsers.parse_key("vlbi_master_file", file_vars=file_vars)
