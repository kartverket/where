"""Makes the radio source coordinates from VASCC available

Description:
------------




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# External library imports
import numpy as np

# Where imports
from where.apriori import crf
from where.lib import files
from where.lib import log
from where.lib import plugins
from where import parsers
from where.lib import rotation
from where.ext import sofa_wrapper as sofa
from where.lib.time import Time
from where.lib.unit import unit


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
        return np.array([source_info["ra"], source_info["dec"]])
