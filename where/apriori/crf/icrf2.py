"""A module to provide information from radio sources defined in ICRF2

Description:
------------

Reads source positions ICRF2 files. Source positions are considered constant in ICRF2. 

References:
-----------
:cite:'icrf2009'

"""

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.data.direction import Direction
from where.apriori import crf
from where import parsers


@plugins.register
class Icrf2(crf.CrfFactory):
    """A class to provide information from radio sources defined in ICRF2
    """

    def _read_data(self):
        """Read data needed by this Celestial Reference Frame for calculating positions of sources

        Returns:
            Dict:  Dictionary containing data about each source defined in this reference frame.
        """
        data = parsers.parse_key(file_key="icrf2_non_vcs").as_dict()
        data.update(parsers.parse_key(file_key="icrf2_vcs_only").as_dict())

        return data

    def _calculate_pos_crs(self, source):
        """Calculate position for a source

        Args:
            source (String):    Key saying which source to calculate position for.

        Returns:
            Array:  Positions, one 2-vector
        """
        source_info = self.data[source]
        vector = Direction(ra=source_info["ra"], dec=source_info["dec"], time=self.time)
        return np.squeeze(vector)
