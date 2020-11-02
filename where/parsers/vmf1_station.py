"""A parser for reading station VMF1 files

Description:
------------

columns:
+ station name (8 characters)
+ modified Julian date
+ hydrostatic coefficient a
+ wet coefficient a
+ hydrostatic zenith delay
+ wet zenith delay
+ mean temperature of the atmosphere above the site in Kelvin
+ pressure at the site in hPa
+ temperature at the site in Celsius
+ water vapour pressure at the site in hPa
+ orthometric height of the station (using geoid EGM96)

Example:
HAYSTACK  43874.00  0.00123547  0.00061937  2.3073  0.1400  270.5  1013.02     3.49   7.26  145.6

References:
-----------

"""
#External library import
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class Vmf1StationParser(LineParser):
    """A parser for reading station dependent VMF1 files
    """

    def setup_parser(self):
        return dict(
            names="station, mjd, ah, aw, zh, zw, mean_atm_temp, pressure, temp, wvp, ortho_height ",
            dtype=("U8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8", "f8"),
            autostrip=True,
        )

    def structure_data(self):
        fieldnames = list(self._array.dtype.names)
        key_field = fieldnames.pop(0)
        keys = np.unique(self._array[key_field])
        for key in keys:
            idx = self._array[key_field] == key
            self.data.setdefault(key, {})
            for field in fieldnames:
                self.data[key][field] = list(self._array[field][idx])

