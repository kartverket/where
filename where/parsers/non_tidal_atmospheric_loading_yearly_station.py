"""A parser for reading 3D displacements due to non tidal atmospheric loading

Description:
------------

Supports the format from https://vmf.geo.tuwien.ac.at/APL_products/VLBI/yearly/
"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class NonTidalAtmosphericLoadingDisplacementsYearlyStationParser(LineParser):
    """A parser for reading 3D displacements due to non tidal atmospheric loading
    """

    def setup_parser(self):
        return dict(
            names="station, mjd, up, east, north",
            dtype=("U8", "f8", "f8", "f8", "f8"),
            autostrip=True,
            comments="!"
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
