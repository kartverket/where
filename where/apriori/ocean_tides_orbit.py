"""Get coefficients for ocean tides

Description:

Reads ocean tides coefficients from file.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Where imports
from where.lib import config
from where import parsers
from where.lib import plugins


@plugins.register
def get_ocean_tides():
    """Get ocean tidal loading coefficients

    Reads ocean tidal coefficients from file using OceanTidesFes2004Parser for satellitte
    displacements.

    Returns:
        A dictionary with information about ocean tidal coefficients.
    """
    model = config.tech.orbit_ocean_tides.str
    file_key = "ocean_tides_{}".format(model) if model else "ocean_tides"

    return parsers.parse(file_key=file_key)
