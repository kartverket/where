"""Get coefficients for ocean tidal loading

Description:

Reads ocean tidal loading coefficients from file.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Where imports
from where.lib import plugins
from where import parsers
from where.lib import config


@plugins.register
def get_ocean_tides():
    """Get ocean tidal loading coefficients

    Reads ocean tidal loading from file using OceanTidesParser for station
    displacements.

    Returns:
        A dictionary with information about ocean tidal loading coefficients.
    """
    ocean_tides_model = config.tech.ocean_tides.str
    file_key = "ocean_tides_{}".format(ocean_tides_model) if ocean_tides_model else "ocean_tides"

    return parsers.parse_key(file_key=file_key).as_dict()
