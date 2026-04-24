""" Test :mod:`where.lib.gnss`.

The test set up compares results from Where against gLAB solution for satellite G07 and epoch 2016-03-01 00:00:00.0,
whereby precise orbits are used.



"""

import datetime
import unittest

import numpy as np

from where.lib import gnss
from where import data


class TestGnss(unittest.TestCase):
    def setUp(self):
        """

        """
        # Definition of date
        year = 2016
        month = 3
        day = 1
        hour = 0
        minute = 0
        second = 0

        rundate = datetime.date(year, month, day)
        time = [
            (
                "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
                "".format(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
            )
        ]

        self.dset = data.Dataset(rundate, tech=None, stage=None, dataset_name="test_where", dataset_id=0, empty=True)

        # Add 'satellite' field to Dataset
        self.dset.num_obs = 1
        self.dset.add_text("satellite", val=list(["G07"]))

        # Add 'time' field to Dataset
        self.dset.add_time("time", val=list([2457448.5]), scale="gps", format="jd")

        # Add 'sat_posvel' field to Dataset
        itrs_pos = np.array([-5499488.571916, -14614558.207892, 21564722.457044])  # [m]
        itrs_vel = np.array([2162.8060, -1636.5180, -599.2837])  # [m/s]
        itrs = np.hstack((itrs_pos, itrs_vel)).reshape(1, 6)
        self.dset.add_posvel("sat_posvel", time="time", itrs=itrs)

        # Add 'site_pos' field to Dataset
        itrs = np.array([[3275753.7411, 321110.9298, 5445042.0022]])  # [m]

        self.dset.add_position("site_pos", time="time", itrs=itrs)

    def test_get_flight_time(self):
        flight_time = gnss.get_flight_time(self.dset)
        expected_flight_time = np.array([0.07893])
        np.testing.assert_allclose(flight_time, expected_flight_time, rtol=0, atol=1e-6)
        # print('OUTPUT:\n flight_time = {:f} [s]\n'.format(flight_time[0]))

    def test_get_satellite_phase_center_offset(self):
        """
        Test against results from gLAB function satellitePhaseCenterCorrection3D() in model.c for satellite G07 and
        epoch 2016-03-01 00:00:00.0 with the ANTEX PCO's in North, East and Up: -0.000400 0.005000 0.822400 m.
        """
        # TODO: move test to apriori/tests
        pco_itrs = gnss.get_satellite_phase_center_offset(self.dset)
        expected_pco_itrs = np.array([[0.169682, 0.455598, -0.663329]])
        np.testing.assert_allclose(pco_itrs, expected_pco_itrs, rtol=0, atol=1e-4)
        # print('OUTPUT:\n pco_itrs = {:f} {:f} {:f} [m]\n'.format(pco_itrs[0][0], pco_itrs[0][1], pco_itrs[0][2]))


if __name__ == "__main__":
    unittest.main()
