"""
Description:

    Reads station-dependent information from file
    which contains informaton on which data sources
    are reliable or not.

References:

    http://ilrs.dgfi.tum.de/fileadmin/data_handling/ILRS_Data_Handling_File.snx



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# Where imports
from where import parsers
from where.lib import log
from where.lib import plugins


@plugins.register
def get_handling_file(time):
    """Read station-dependent info
    """
    handling = parsers.parse(file_key="slr_handling_file", time=time)
    for station in handling:
        if "V" in handling[station]:
            log.warn("New station {}, unreliable coordinates according to ILRS Handling file", station)
        if "E" in handling[station]:
            log.warn("Should esimate bias for station {}, according to ILRS Handling file", station)
        if "R" in handling[station]:
            log.warn("Should apply bias for station {}, according to ILRS Handling file", station)
        if "U" in handling[station]:
            log.warn("Should estimate time bias for station {}, according to ILRS Handling file", station)
        if "X" in handling[station]:
            log.warn("Should delete data for station {}, according to ILRS Handling file", station)
    return handling
