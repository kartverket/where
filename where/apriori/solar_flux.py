"""Get a priori solar flux data

Description:
------------

Reads time-dependent solar flux from file.

References:
-----------

http://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-radio/noontime-flux/penticton/penticton_observed/tables/

https://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-radio/noontime-flux/penticton/documentation/dataset-description_penticton.pdf

TODO: All these files from the link above are merged to F10.7CM (leftovers from old GEOSAT). Do something about this?



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""
from datetime import timedelta
from scipy import interpolate

# Where imports
from where import parsers
from where.lib import config
from where.lib import plugins


@plugins.register
def get_solar_flux(rundate):
    """Read time-dependent solar flux from file
    """
    flux = parsers.parse(file_key="solar_flux")
    arc_length = config.tech.arc_length.int

    date_to_add = rundate - timedelta(days=5)
    time_list, flux_list = list(), list()

    while True:
        if date_to_add in flux:
            # Note that 72000 sec corresponds to local noontime in Penticton
            time_list.append((date_to_add - rundate).total_seconds() + 72000)
            flux_list.append(flux[date_to_add])
            if date_to_add > rundate + timedelta(days=arc_length):
                break
        date_to_add += timedelta(days=1)

    return interpolate.interp1d(time_list, flux_list)
