"""Get apriori data for eccentricities

Description:

Reads data for eccentricity files




"""
# Where imports
from where import parsers
from where.lib import plugins


@plugins.register
def get_eccentricity(rundate):
    """Get Eccentricities for a given date

    Args:
        rundate: The run date of the data analysis.

    Returns:
        A dictionary of eccentricities.
    """

    return parsers.parse(file_key="eccentricity", rundate=rundate)
