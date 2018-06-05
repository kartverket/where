"""Calculate the station displacements due to ocean tidal loading

Description:
------------

Correct for ocean tides. This model uses the IERS Software Collection HARDISP-function to calculate the local
displacement of a station due to ocean tidal loading. The displacement is returned in a geocentric celestial system.


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] Lyard, F., Lefevre F., Letellier T. and Francis, O,
       Modelling the global ocean tides: modern insights from FES2004,
       Ocean Dynamics (2006) 56: 394-415.
       http://maia.usno.navy.mil/conv2010/chapter6/add_info/FES2004.pdf

.. [3] Biancale, R., Ocean and atmospheric tides standards,
       Presentation at IERS Workshop on Conventions, September 2007.
       http://maia.usno.navy.mil/conv2010/chapter6/add_info/BIPM_ocean_tides_2007.pdf




$Revision: 15044 $
$Date: 2018-05-08 23:23:38 +0200 (Tue, 08 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.ext import iers_2010 as iers
from where.lib import config
from where.lib import log
from where.lib import plugins


@plugins.register
def ocean_tides(dset):
    """Calculate ocean tide corrections at both stations

    Ocean tide corrections are returned in meters in the Geocentric Celestial Reference System for each observation. A
    Numpy array with 6 columns is returned, the first three columns are \f$ x, y, z \f$ for station 1, while the last
    three columns are \f$ x, y, z \f$ for station 2.

    Args:
        rundate:  The model run date.
        tech:     The technique.
        dset:     A Dataset containing model data.

    Returns:
        Numpy array with ocean tide corrections in meters.
    """
    # Ocean Tides Coefficients
    otc = apriori.get("ocean_tides")

    amplitudes = {s: otc[s]["amplitudes"] for s in dset.unique("site_id") if s in otc}
    phases = {s: otc[s]["phases"] for s in dset.unique("site_id") if s in otc}
    # Warn about missing Ocean Tides Coefficients
    for site_id in set(dset.unique("site_id")) - set(amplitudes.keys()):
        log.error("Missing ocean loading coefficients for site id '{}'. Correction set to zero.", site_id)

    data_out = list()
    for _ in dset.for_each("station"):
        data_out.append(ocean_tides_station(dset, amplitudes, phases))

    return np.hstack(data_out)


def ocean_tides_station(dset, amplitudes, phases):
    """Calculate the ocean tide corrections for a station

    Ocean tide corrections are returned in meters in the Geocentric Celestial Reference System for each observation.

    Args:
        dset:        A Dataset containing model data.

    Returns:
        Numpy array with ocean tide corrections in meters.
    """
    denu = np.zeros((dset.num_obs, 3))
    use_cmc = config.tech.ocean_tides_cmc.bool

    # Calculate correction
    for obs, (site_id,) in enumerate(dset.values("site_id")):
        if site_id not in amplitudes:
            continue

        epoch = [float(t) for t in dset.time.utc[obs].yday.split(":")]
        dup, dsouth, dwest = iers.hardisp(epoch, amplitudes[site_id], phases[site_id], 1, 1.0)

        # Correction in topocentric (east, north, up) coordinates
        denu[obs] = np.array([-dwest[0], -dsouth[0], dup[0]])

    dxyz = dset.site_pos.convert_enu_to_itrs(denu)

    # Center of mass corrections
    if use_cmc:
        coeff_cmc = apriori.get("ocean_tides_cmc")
        in_phase = coeff_cmc["in_phase"]
        cross_phase = coeff_cmc["cross_phase"]
        for obs, time in enumerate(dset.time.utc):
            year, doy = time.datetime.year, float(time.datetime.strftime("%j")) + time.mjd_frac
            angle = iers.arg2(year, doy)[:, None]
            dxyz[obs] += np.sum(in_phase * np.cos(angle) + cross_phase * np.sin(angle), axis=0)

    return dset.site_pos.convert_itrs_to_gcrs(dxyz)
