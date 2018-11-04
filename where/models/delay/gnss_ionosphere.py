"""Calculate the delay caused by the ionospheric refraction

Description:
------------
"""
# External library imports
import numpy as np

# Midgard imports
from midgard.ionosphere import klobuchar

# Where imports
from where import apriori
from where.lib import log
from where.lib import plugins
#from glpy import gnss_iono_models


@plugins.register
def gnss_ionosphere(dset):
    """Determine ionosphere correction

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray: Ionosphere correction in meter
    """

    # TODO_Mohammed: NeQuick and Klobuchar model is available via glpy.gnss_iono_models. See also klobuchar.py and
    #               nequickg.py example under ./glpy/glpy/examples.
    logger = lambda _: None
        
    orbit = apriori.get(
            "orbit",
            rundate=dset.rundate,
            time=dset.time,
            satellite=tuple(dset.satellite),
            system=tuple(dset.system),
            station=dset.vars["session"],
        )
    rundate = dset.rundate.strftime("%Y-%m-%d")
    iono_para = orbit.dset_raw.meta[rundate]["iono_para"]
    iono_alpha = iono_para["GPSA"]["para"]
    iono_beta = iono_para["GPSB"]["para"]
    rec_pos = dset.site_pos.llh.mean(axis=0)

    corrections = np.zeros(dset.num_obs) 
    for idx, (time, az, el) in enumerate(dset.values("time", "site_pos.azimuth", "site_pos.elevation")):
        try:
            delay, _ = klobuchar.klobuchar(time.gpssec, iono_alpha + iono_beta, rec_pos, az, el, logger=logger)
        except ValueError as err:
            log.warn(f"Problem with ionosphere model: {err}")
        corrections[idx] = delay
    import IPython; IPython.embed()
    
        
    return corrections
