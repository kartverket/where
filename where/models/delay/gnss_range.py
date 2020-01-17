"""Calculate the station-satellite distance

Description:
------------

This model computes the distance between the station and the satellite.

"""
# Midgard imports
from midgard.dev import plugins


@plugins.register
def gnss_range(dset: "Dataset"):
    """Calculate distance between station and satellite in GCRS

    Args:
        rundate:    The model run date.
        tech:       Name of technique.
        dset:       A dataset containing the data.

    Returns:
        table of corrections for each observation
    """
    # If the range between the satellite and the receiver is determined in the ITRS, then it has to be taken into
    # account that the Earth is rotating during the time of flight of the satellite signal. This is not necessary in
    # the GCRS.
    correction = (dset.sat_posvel.gcrs.pos - dset.site_pos.gcrs.pos).length

    # if "gnss_earth_rotation" in dset.fields:
    #    dset.sat_posvel[:] += dset.gnss_earth_rotation
    #    correction_itrs = (dset.sat_posvel.trs.pos - dset.site_pos.trs.pos).length
    #    print("MURKS: ", correction_itrs - correction)
    #    correction = correction_itrs

    return correction
