"""A parser for reading 3D displacements due to non tidal atmospheric loading

Description:
------------

Supports the format from http://vmf.geo.tuwien.ac.at/APL_products/GRID
"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_line import LineParser


@plugins.register
class NonTidalAtmosphericLoadingDisplacementsParser(LineParser):
    """A parser for reading 3D displacements due to non tidal atmospheric loading
    """

    def setup_parser(self):
        return dict(names=["lat", "lon", "up", "east", "north"], comments="!")

    def structure_data(self):
        # Convert displacements from array to matrix
        lat = np.unique(self._array["lat"])
        lon = np.unique(self._array["lon"])
        num_value = (len(lat), len(lon))

        # Latitude is given from 90 degrees to -90 degrees, reverse order
        # The latitude angles themselves are already sorted by np.unique
        up = self._array["up"].reshape(num_value)[::-1, :]
        east = self._array["east"].reshape(num_value)[::-1, :]
        north = self._array["north"].reshape(num_value)[::-1, :]

        # Shift longitude -180 degrees
        lon = (lon + 180) % 360 - 180
        idx = lon.argsort()
        lon = lon[idx]
        up = up[:, idx]
        east = east[:, idx]
        north = north[:, idx]

        lat = np.radians(lat)
        lon = np.radians(lon)

        self.data = dict(lat=lat, lon=lon, up=up, east=east, north=north)
