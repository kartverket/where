""" Test :mod:`where.models.delay.gnss_satellite_phase_center_offset`.


DATASET:

Tests are compared against results from gLAB 3.0.0 version for day 2016-03-01 and station STAS.

itrs_pos = (-5499488.548910, -14614558.223893, 21564722.455041)  # satellite position vector in ITRS in [m]
itrs_vel = (2162.8060, -1636.5180, -599.2837)                    # satellite velocity vector in ITRS in [m/s]
site_pos = (3275753.7411, 321110.9298, 5445042.0022)             # station position in ITRS in [m]

Satellite phase center correction returned by routine satellitePhaseCenterCorrection() in model.c:
    pc = -0.802373 m


$Revision: 15940 $
$Date: 2018-09-11 16:32:58 +0200 (Tue, 11 Sep 2018) $
$LastChangedBy: hjegei $
"""

# Standard library imports
from datetime import date
import unittest

# External library imports
import numpy as np

# Where imports
from where import data
from where.models.delay import gnss_satellite_phase_center_offset
from where.lib import log


class TestPositionTable(unittest.TestCase):
    def setUp(self):

        # Print debugging messages
        # log.init('DEBUG')

        # Definition of date
        year = 2016
        month = 3
        day = 1

        rundate = date(year, month, day)
        time = [("{year}-{month:02d}-{day:02d}" "".format(year=year, month=month, day=day))]

        self.dset = data.Dataset(rundate, tech=None, stage=None, dataset_name="test_where", dataset_id=0, empty=True)

        # Add 'satellite' field to Dataset
        self.dset.num_obs = 1
        self.dset.add_text("satellite", val=list(["G07"]))

        # Add 'time' field to Dataset
        self.dset.add_time("time", val=list([2457448.5]), scale="gps", format="jd")

        # Add 'sat_posvel' field to Dataset
        itrs_pos = np.array([-5499488.548910, -14614558.223893, 21564722.455041])  # [m]
        itrs_vel = np.array([2162.8060, -1636.5180, -599.2837])  # [m/s]
        itrs = np.hstack((itrs_pos, itrs_vel)).reshape(1, 6)
        self.dset.add_posvel("sat_posvel", time="time", itrs=itrs)

        # Add 'site_pos' field to Dataset
        itrs = np.array([[3275753.7411, 321110.9298, 5445042.0022]])  # [m]
        self.dset.add_position("site_pos", time="time", itrs=itrs)

    def test_gnss_satellite_phase_center_offset(self):
        correction = gnss_satellite_phase_center_offset.gnss_satellite_phase_center_offset(self.dset)

        expected_correction = np.array([-0.802373])
        np.testing.assert_allclose(correction, expected_correction, rtol=0, atol=1e-4)
        # print('OUTPUT:\n correction = {:f} [m]'.format(-correction[0]))


if __name__ == "__main__":
    unittest.main()
