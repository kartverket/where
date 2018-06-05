"""Get apriori data for eccentricities

Description:

Reads data for eccentricity files



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

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
