"""A parser for reading data from ITRF files in SSC format

Description:
------------

Reads station positions and velocities from ITRF files in SSC format. The velocity model is a simple linear offset
based on the reference epoch.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# External library imports
import numpy as np

# Where imports
from where.apriori import trf
from where.lib import plugins
from where.lib.time import Time
from where.lib.unit import unit
from where import parsers


@plugins.register
class Vascc(trf.TrfFactory):
    """A class for representing apriori station positions and velocities from SSC
    """

    def _read_data(self):
        return parsers.parse_key(file_key="vascc_trf").as_dict()

    #
    # Velocity model
    #
    def _calculate_pos_itrs(self, site):
        """Calculate positions for the given time epochs

        The positions are calculated as simple linear offsets based on the reference epoch. Makes sure to pick out the
        correct time interval to use.

        Args:
            key:    Key saying which site to calculate position for, type might depend on the Trf.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        station_info = self.data[site]
        ref_epoch = Time(float(station_info["ref_epoch"]), format="decimalyear", scale="utc")
        pos = np.full((self.time.size, 3), fill_value=np.nan)

        ref_pos = np.array(station_info["pos"])
        ref_vel = np.array(station_info["vel"])
        interval_years = (self.time - ref_epoch).jd * unit.day2julian_years
        pos[:, :] = ref_pos + interval_years[:, None] * ref_vel[None, :]

        return pos
