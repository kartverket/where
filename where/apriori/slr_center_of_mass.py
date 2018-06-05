"""Description:

    Reads station-dependent information from file.  This is needed in order to compute center of mass corrections for
    the satellites lageos-1, lageos-2, etalon-1, etalon-2 and ajisai.

References:

    http://ilrs.dgfi.tum.de


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

# Where imports
from where.lib import plugins
from where import parsers


@plugins.register
def get_center_of_mass(sat_name):
    """Read station-dependent center of mass corrections from file

    Args:
        sat_name (String):  Name of satellite.

    Returns:
        Dict:  Center of mass corrections per station for the given satellite.
    """
    return parsers.parse("slr_com", sat_name=sat_name)
