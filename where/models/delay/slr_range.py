"""Calculate the station-satellite distance

Description:
------------

This is the main SLR model, computing the distance between the station and the satellite, taking into account the
gravity field of the Earth.



$Revision: 14978 $
$Date: 2018-04-30 19:01:11 +0200 (Mon, 30 Apr 2018) $
$LastChangedBy: hjegei $

"""
# Where imports
from where.lib import constant
from where.lib import plugins


@plugins.register
def slr_range(dset):
    """Calculate the distance between station and satellite

    Integrate differential equation of motion of the satellite and differential equation of the state transition matrix
    of the satellite.

    Args:
        dset:       A dataset containing the data

    Returns:
        Numpy array: Distance for each observation in meters
    """
    return (dset.up_leg + dset.down_leg) / 2 * constant.c
