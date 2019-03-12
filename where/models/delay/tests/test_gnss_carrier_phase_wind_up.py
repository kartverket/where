""" Test :mod:`where.models.delay.gnss_carrier_phase_wind_up`.

The test set up compares results from Where against gLAB solution for satellite G07 and epoch 2016-03-01 00:00:00.0,
whereby precise orbits are used.

$Revision: 15940 $
$Date: 2018-09-11 16:32:58 +0200 (ti., 11 sep. 2018) $
$LastChangedBy: hjegei $
"""

# Standard library imports
from datetime import datetime
import unittest

# External library imports
import numpy as np

# Where imports
from where import data
from where.models.delay import gnss_carrier_phase_wind_up


class TestPositionTable(unittest.TestCase):
    def setUp(self):

        # Definition of date
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

    def test_gnss_carrier_phase_wind_up(self):
        """
        Result: -6.00897
        """
        correction = gnss_carrier_phase_wind_up.gnss_carrier_phase_wind_up(self.dset)

        # expected_correction = np.array([[-6.00897]])

        # np.testing.assert_allclose(correction, expected_correction, rtol=0, atol=1e-5)

        print("OUTPUT:\n correction = {:f} [m]".format(correction[0]))


#    def test_get_yaw_coord_sys(self):
#        """
#        Test is done against following results of gLAB_3.0.0 software package (given by model.c routines
#        fillSatelliteOrientation() and getSatelliteOrientation()):

#        INPUT:
#            Station:                              STAS (Stavanger)
#            Modified Julian Day:                  57448.0  (2013-03-01)
#            GPS PRN:                              7
#            Sun position vector in ECEF:          (-148102175.801620, -7980320.884334, -19534142.160482) [m]
#            Satellite position vector in ECEF:    (-5499488.548910, -14614558.223893, 21564722.455041) [m]
#
#        OUTPUT:
#            unit_x vector: (-0.971554, 0.017062, -0.236205) [m]
#            unit_y vector: (-0.115836, 0.835705,  0.536822) [m]
#            unit_z vector: ( 0.206557, 0.548913, -0.809956) [m]

#        Note: gLAB does not use JPL ephemeris for determination of geocentric Sun position vector, so differences should
#              be expected.
#        """
#        unit_x, unit_y, unit_z = gnss_carrier_phase_wind_up.get_yaw_coord_sys(self.dset.time,
#                                                                              self.dset.sat_posvel.itrs_pos)

#        expected_unit_x = np.array([[ -0.971554, 0.017062, -0.236205]])
#        expected_unit_y = np.array([[ -0.115836, 0.835705, 0.536822]])
#        expected_unit_z = np.array([[ 0.206557, 0.548913, -0.809956]])
#
#        np.testing.assert_allclose(unit_x, expected_unit_x, rtol=0, atol=1e-5)
#        np.testing.assert_allclose(unit_y, expected_unit_y, rtol=0, atol=1e-5)
#        np.testing.assert_allclose(unit_z, expected_unit_z, rtol=0, atol=1e-5)

#        #print('OUTPUT:\n unit_x = {:f} {:f} {:f} [m]'.format(unit_x[0,0], unit_x[0,1], unit_x[0,2]))
#        #print(' unit_y = {:f} {:f} {:f} [m]'.format(unit_y[0,0], unit_y[0,1], unit_y[0,2]))
#        #print(' unit_z = {:f} {:f} {:f} [m]'.format(unit_z[0,0], unit_z[0,1], unit_z[0,2]))


if __name__ == "__main__":
    unittest.main()
