"""Makes the radio source coordinates from VASCC available

Description:
------------





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
class Vascc(crf.CrfFactory):
    """A class to provide information from radio sources defined in ICRF2
    """

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.format> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        return parsers.parse_key(file_key="vascc_crf").as_dict()

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
