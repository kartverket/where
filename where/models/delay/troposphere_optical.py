"""Calculate total zenith delay

Description:
------------

This model determines the total zenith delay following Mendes and Pavlis [1].

References:
-----------

[1] Mendes, V.B. and E.C. Pavlis, 2004,
    "High-accuracy zenith delay prediction at optical wavelengths,"
    Geophysical Res. Lett., 31, L14602, doi:10.1029/2004GL020308, 2004

[2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010)



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""
# External library imports
import numpy as np

# Where imports
from where.ext import iers_2010 as iers
from where.lib import plugins
from where.lib.unit import unit


@plugins.register
def pavlis_mendes(dset):
    """Calculate zenith delay for all observations

    Args:
        dset:     Dataset containing the data

    Returns:
        Numpy array:  Correction in meters for each observation
    """
    output = np.zeros(dset.num_obs)

    # Compute WVP = Water Vapour Pressure from temperature and humidity:
    #   https://en.wikipedia.org/wiki/Vapour_pressure_of_water
    wvp = np.exp(20.386 - 5132 / dset.temperature) * unit.mmHg2hPa * dset.humidity * unit.percent2unit

    for obs, ((lat, _, height), pressure, wavelength, temperature, elevation) in enumerate(
        dset.values("site_pos.llh", "pressure", "wavelength", "temperature", "site_pos.elevation")
    ):
        # Compute total zenith delay:
        output[obs] = iers.fculzd_hpa(
            np.degrees(lat), height, pressure, wvp[obs], wavelength * unit.nanometer2micrometer
        )[
            0
        ]
        # Mapping function:
        output[obs] *= iers.fcul_a(np.degrees(lat), height, temperature, np.degrees(elevation))

    return output
