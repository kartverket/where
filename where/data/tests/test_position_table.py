""" Test :mod:`where.data.position_table`.

-------


"""

# Standard library imports
from datetime import datetime
import unittest

# External library imports
import numpy as np

# Where imports
from where import data


class TestPositionTable(unittest.TestCase):
    def setUp(self):

        # Definition of dummy date
        year = 2016
        month = 3
        day = 1
        hour = 0
        minute = 0
        second = 0

        rundate = datetime(year, month, day, hour, minute, second)
        time = [
            (
                "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
                "".format(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
            )
        ]

        self.dset = data.Dataset(rundate, tech=None, stage=None, dataset_name="test_where", dataset_id=0, empty=True)

        # Add 'site_pos' field to Dataset
        itrs = np.array([[3771793.968, 140253.342, 5124304.349]])  # [m]
        self.dset.num_obs = 1
        self.dset.add_time("time", val=time, scale="utc")
        self.dset.add_position("site_pos", time="time", itrs=itrs)

    def test_llh(self):
        """
        Based on example in Section 2.2.1 of :cite:`ogp2011` (Note: The example uses WGS84 reference ellipsoid instead
        of GRS80.)

        Result:
           latitude  = 53 [deg] 48' 33.820'' = 53.80939444444444 [deg] = 0.9391511015599004 [rad]
           longitude =  2 [deg] 07' 46.380'' = 2.12955 [deg] = 0.037167659085845246 [rad]
           height = 73.0 m
        """
        llh = self.dset.site_pos.llh
        expected_llh = np.array([[0.93915110, 0.03716766, 73.0]])  # station position in geodetic coordinates (
        # latitude [rad], longitude [rad], height [m])

        np.testing.assert_allclose(llh, expected_llh, rtol=0, atol=1e-5)

        # print('OUTPUT:\n llh = {:f} [rad] {:f} [rad] {:f} [m]\n'
        #       ''.format(llh[0,0], llh[0,1], llh[0,2]))

    def test_convert_itrs_to_enu(self):
        """
        Geocentric coordinates from argument `itrs` are related to the geocenter, which is also valid for the
        "topocentric" (in this case "geocentric") coordinates East, North and Up.
        """
        itrs = np.array([[-329.4395, 167.8142, 9.2296]])
        enu = self.dset.site_pos.convert_itrs_to_enu(itrs)
        expected_enu = np.array([[179.94, 266.11, -183.26]])  # station position in East, North and Up in [m] related
        # to geocenter

        np.testing.assert_allclose(enu, expected_enu, rtol=0, atol=1e-4)

        # print('OUTPUT:\n enu = {:f} [m] {:f} [m] {:f} [m]\n'
        #       ''.format(enu[0,0], enu[0,1], enu[0,2]))

    def test_convert_enu_to_itrs(self):
        """
        Geocentric coordinates from argument `itrs` are related to the geocenter, which is also valid for the
        "topocentric" (in this case "geocentric") coordinates East, North and Up.
        """
        enu = np.array([[179.94, 266.11, -183.26]])
        itrs = self.dset.site_pos.convert_enu_to_itrs(enu)
        expected_itrs = np.array([[-329.4395, 167.8142, 9.2296]])  # position in X, Y and Z in [m]
        np.testing.assert_allclose(itrs, expected_itrs, rtol=0, atol=1e-4)

        # print('OUTPUT:\n itrs = {:f} [m] {:f} [m] {:f} [m]\n'
        #       ''.format(itrs[0,0], itrs[0,1], itrs[0,2]))


if __name__ == "__main__":
    unittest.main()
