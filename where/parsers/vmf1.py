"""A parser for reading VMF1 files

Description:
------------

Reads the gridded coefficient 'a' for the wet and hydrostatic delay model.
Reads the gridded zenith wet and hydrostatic delay from the model.
Reads the gridded ellipsoidal heights.

The first line in each VMF1 gridded file is a header showing the values in degrees for north, south, west, east,
spacing north-south, spacing west-east. The rest of the file contains the parameters in latitude bands going from north
to south (90 to -90 degrees) in 2.0 degree steps, and from west to east within each band (0 to 360 degrees) in 2.5
degree steps.

References:
-----------
http://ggosatm.hg.tuwien.ac.at/DELAY/readme.txt

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser


@plugins.register
class Vmf1Parser(Parser):
    """A parser for reading VMF1 files
    """

    def read_data(self):
        data = np.fromfile(str(self.file_path), sep=" ")
        num_value = (int(abs(data[1] - data[0]) / data[4] + 1), int(abs(data[3] - data[2]) / data[5] + 1))
        lat = np.linspace(data[0], data[1], num_value[0])
        lon = np.linspace(data[2], data[3], num_value[1])
        values = np.array(data[6:]).reshape(num_value)

        # Latitude is given from 90 degrees to -90 degrees, reverse order
        lat = lat[::-1]
        values = values[::-1, :]

        # Strip the last longitude to avoid double entry for 0 degrees
        lon = lon[:-1]
        values = values[:, :-1]

        # Shift longitude -180 degrees
        lon = (lon + 180) % 360 - 180
        idx = lon.argsort()
        lon = lon[idx]
        values = values[:, idx]

        lat = np.radians(lat)
        lon = np.radians(lon)
        self.data = dict(lat=lat, lon=lon, values=values)
