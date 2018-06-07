""" Test :mod:`where.data.posvel_table`.

Description:
------------

DATASET 1:

Several tests are based on exercise 2.3 in :cite:`montenbruck2012` with following example:

r = (10000000.0, 40000000.0, -5000000.0)  # satellite position vector in [m]
v = (-1500.0, 1000.0, -100.0)             # satellite velocity vector in [m/s]

a     = 25015181.0181 [m]   # Semi-major axis
e     = 0.7079772           # Eccentricity
i     =   6.971 [deg]       # Inclination
Omega = 173.290 [deg]       # Right ascension of ascending node
omega =  91.553 [deg]       # Argument of perigee
M     = 144.225 [deg]       # Mean anomaly



DATASET 2:

Other tests are compared against results from gLAB 3.0.0 version  for day 2016-03-01 and satellite G07
(epoch: 57448.0 MJD).

r = (-5499233.574326, -14614718.575397, 21564674.490672)  # satellite position vector in [m]
v = (3228.525601, -2037.522778, -599.319641)              # satellite velocity vector in [m/s]

Rotation matrix returned by routine getSatelliteOrientationACR() in model.c:

    [ ea ]   [  0.833664, -0.532072, -0.148000 ]   # Rotation matrix from geocentric XYZ coordinates to local orbital
R = [ ec ] = [  0.512194,  0.644661,  0.567512 ]   # reference system (along-track, cross-track and radial)
    [ er ]   [ -0.206547, -0.548919,  0.809954]

DATASET 3:

Other tests are compared against results from gLAB 3.0.0 version for day 2016-03-01 and station STAS (satellite: G07,
epoch: 57448.0 MJD) by using PPP solution (with precise orbits).

r = (-5499488.571916, -14614558.207892, 21564722.457044)  # satellite position vector in [m]
v = (0, 0, 0)  # satellite velocity vector in [m/s]

Rotation matrix returned by routine getSatelliteOrientation() in model.c:

    [ ex ]   [ -0.971554, 0.017062, -0.236205 ]   # Rotation matrix from geocentric XYZ coordinates to yaw-steering
R = [ ey ] = [ -0.115836, 0.835705,  0.536822 ]   # reference system
    [ ez ]   [  0.206557, 0.548913, -0.809956 ]

and the belonging unit vector pointing from the satellite to Sun is r_sun = (-0.959930, 0.044658, -0.276657) m


"""

# Standard library imports
from datetime import datetime
import unittest

# External library imports
import numpy as np

# Where imports
from where import data
from where.lib import gnss
from where.lib import mathp


class TestPosVelTable(unittest.TestCase):

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

        # Define 1st Dataset
        # ------------------
        self.dset = data.Dataset(rundate, tech=None, stage=None, dataset_name="test_where", dataset_id=0, empty=True)

        # Add 'sat_posvel' field to Dataset
        itrs_pos = np.array([10000000.0, 40000000.0, -5000000.0])  # [m]
        itrs_vel = np.array([-1500.0, 1000.0, -100.0])  # [m/s]
        itrs = np.hstack((itrs_pos, itrs_vel)).reshape(1, 6)
        self.dset.num_obs = 1
        self.dset.add_time("time", val=time, scale="gps")
        self.dset.add_posvel("sat_posvel", time="time", itrs=itrs)

        # Add 'kepler' array
        self.kepler = np.array(
            [
                [
                    25015181.018465,  # a [m]
                    0.707977170873199,  # e
                    0.121662175957290,  # i [rad]
                    3.024483909022929,  # Omega [rad]
                    1.597899323919624,  # omega [rad]
                    2.772570719534964,
                ]
            ]
        )  # E [rad]

        # Define 2nd Dataset
        # ------------------
        self.dset2 = data.Dataset(rundate, tech=None, stage=None, dataset_name="test_where2", dataset_id=0, empty=True)

        # Add 'sat_posvel' field to Dataset
        itrs_pos = np.array([-5499233.574326, -14614718.575397, 21564674.490672])  # [m]
        itrs_vel = np.array([3228.525601, -2037.522778, -599.319641])  # [m/s]
        itrs = np.hstack((itrs_pos, itrs_vel)).reshape(1, 6)
        self.dset2.num_obs = 1
        self.dset2.add_time("time", val=time, scale="gps")
        self.dset2.add_posvel("sat_posvel", time="time", itrs=itrs)

        # Define 3rd Dataset
        # ------------------
        self.dset3 = data.Dataset(rundate, tech=None, stage=None, dataset_name="test_where3", dataset_id=0, empty=True)

        # Add 'sat_posvel' field to Dataset
        itrs_pos = np.array([-5499488.571916, -14614558.207892, 21564722.457044])  # [m]
        itrs_vel = np.array([0, 0, 0])  # [m/s]
        itrs = np.hstack((itrs_pos, itrs_vel)).reshape(1, 6)
        self.dset3.num_obs = 1
        self.dset3.add_time("time", val=time, scale="tdb")
        self.dset3.add_posvel("sat_posvel", time="time", itrs=itrs)

    def test_convert_kepler_to_itrs(self):
        posvel = self.dset.sat_posvel.convert_kepler_to_itrs(self.kepler)
        expected_posvel = np.array(
            [[10000000.0, 40000000.0, -5000000.0, -1500.0, 1000.0, -100.0]]  # satellite position in [m]
        )  # satellite velocity in [m/s]

        np.testing.assert_allclose(posvel, expected_posvel, rtol=0, atol=1e-6)

        # print('OUTPUT:\n pos = {:f} {:f} {:f} [m]\n vel = {:f} {:f} {:f} [m/s]\n'
        #       ''.format(posvel[0,0], posvel[0,1], posvel[0,2],
        #                 posvel[0,3], posvel[0,4], posvel[0,5]))

    def test_mean_anomaly(self):
        M = self.dset.sat_posvel.mean_anomaly
        expected_M = 2.51720096  # [rad]
        np.testing.assert_allclose(M, expected_M, rtol=0, atol=1e-8)
        # print('OUTPUT: M = {:f} [deg]\n'.format(np.rad2deg(M[0])))

    def test_kepler(self):
        kepler = self.dset.sat_posvel.kepler
        expected_kepler = np.array([[25015181.018464539, 0.70797717, 0.12166218, 3.024483909, 1.59789932, 2.77257072]])
        np.testing.assert_allclose(kepler, expected_kepler, rtol=0, atol=1e-8)

        # print('OUTPUT:\n a = {:30.18f} [m]\n e = {:20.18f}\n i = {:20.18f} [deg]\n Omega = {:20.18f} [deg]\n omega = {:20.18f} [deg]\n E = {:20.18f} [deg]\n'
        #       ''.format(kepler[0,0],
        #                 kepler[0,1],
        #                 kepler[0,2],
        #                 kepler[0,3],
        #                 kepler[0,4],
        #                 kepler[0,5]))

    def test_true_anomaly(self):
        vega = self.dset.sat_posvel.true_anomaly
        expected_vega = 2.98755476  # [rad]
        np.testing.assert_allclose(vega, expected_vega, rtol=0, atol=1e-8)
        # print('OUTPUT: vega = {:f} [deg]\n'.format(np.rad2deg(vega[0])))

    def test_acr2itrs(self):
        acr2itrs = self.dset2.sat_posvel._acr2itrs
        expected_acr2itrs = np.array(
            [[[0.833664, 0.512194, -0.206547], [-0.532072, 0.644661, -0.548919], [-0.148000, 0.567512, 0.809954]]]
        )
        np.testing.assert_allclose(acr2itrs, expected_acr2itrs, rtol=0, atol=1e-5)
        # print('OUTPUT:\n acr2itrs = {}'.format(acr2itrs))

    def test_itrs2acr(self):
        itrs2acr = self.dset2.sat_posvel._itrs2acr
        expected_itrs2acr = np.array(
            [[[0.833664, -0.532072, -0.148000], [0.512194, 0.644661, 0.567512], [-0.206547, -0.548919, 0.809954]]]
        )
        np.testing.assert_allclose(itrs2acr, expected_itrs2acr, rtol=0, atol=1e-5)
        # print('OUTPUT:\n itrs2acr = {}'.format(itrs2acr))

    def test_convert_acr_to_itrs(self):
        acr = np.array([[-0.732580, 0.048866, -0.639682]])  # posDiff variable in compareOrbits() routine in gLAB.c
        itrs = self.dset2.sat_posvel.convert_acr_to_itrs(acr)
        expected_itrs = np.array([[-0.453572, 0.772421, -0.381959]])
        np.testing.assert_allclose(itrs, expected_itrs, rtol=0, atol=1e-5)
        # print('OUTPUT:\n itrs = {}'.format(itrs))

    def test_convert_itrs_to_acr(self):
        itrs = np.array([[-0.453572, 0.772421, -0.381959]])
        acr = self.dset2.sat_posvel.convert_itrs_to_acr(itrs)
        expected_acr = np.array([[-0.732580, 0.048866, -0.639682]])
        np.testing.assert_allclose(acr, expected_acr, rtol=0, atol=1e-5)
        # print('OUTPUT:\n acr = {}'.format(acr))

    def test_itrs2yaw(self):
        """
        Test against results from gLAB function getSatelliteOrientation() in model.c.

        gLAB uses simplified expressions for determination of solar coordinates (based on Section 3.3.2 in
        :cite:`montenbruck2012`), whereas Where uses more sophisticated JPL ephemeris.
        """
        itrs2yaw = self.dset3.sat_posvel._itrs2yaw
        expected_itrs2yaw = np.array(
            [
                [
                    [
                        -0.971554, 0.017062, -0.236205
                    ],  # orientation variable in getSatelliteOrientation() routine in model.c
                    [-0.115836, 0.835705, 0.536822],
                    [0.206557, 0.548913, -0.809956],
                ]
            ]
        )
        np.testing.assert_allclose(itrs2yaw, expected_itrs2yaw, rtol=0, atol=1e-2)
        # print('OUTPUT:\n itrs2yaw = {}'.format(itrs2yaw))

    def test_yaw2itrs(self):
        """
        Test against results from gLAB function getSatelliteOrientation() in model.c.

        gLAB uses simplified expressions for determination of solar coordinates (based on Section 3.3.2 in
        :cite:`montenbruck2012`), whereas Where uses more sophisticated JPL ephemeris.
        """
        yaw2itrs = self.dset3.sat_posvel._yaw2itrs
        expected_yaw2itrs = np.array(
            [
                [
                    [
                        -0.971554, -0.115836, 0.206557
                    ],  # orientation variable in getSatelliteOrientation() routine in model.c
                    [0.017062, 0.835705, 0.548913],
                    [-0.236205, 0.536822, -0.809956],
                ]
            ]
        )
        np.testing.assert_allclose(yaw2itrs, expected_yaw2itrs, rtol=0, atol=1e-2)
        # print('OUTPUT:\n yaw2itrs = {}'.format(yaw2itrs))

    def test_convert_yaw_to_itrs(self):
        """
        Test against results from gLAB function satellitePhaseCenterCorrection3D() in model.c.
        """
        yaw = np.array([[-0.000400, 0.005000, 0.822400]])
        itrs = self.dset3.sat_posvel.convert_yaw_to_itrs(yaw)
        expected_itrs = np.array([[0.169682, 0.455598, -0.663329]])
        np.testing.assert_allclose(itrs, expected_itrs, rtol=0, atol=1e-4)
        # print('OUTPUT test_convert_yaw_to_itrs:\n itrs = {} m'.format(itrs))

    def test_convert_itrs_to_yaw(self):
        """
        Test against results from gLAB function satellitePhaseCenterCorrection3D() in model.c.
        """
        itrs = np.array([[0.16968192802350601, 0.45559755989723694, -0.66332935011041794]])
        yaw = self.dset3.sat_posvel.convert_itrs_to_yaw(itrs)
        expected_yaw = np.array([[-0.000400, 0.005000, 0.822400]])
        np.testing.assert_allclose(yaw, expected_yaw, rtol=0, atol=1e-4)
        # print('OUTPUT test_convert_itrs_to_yaw:\n yaw = {} m \n\n\n'.format(yaw))

    def test_itrs_pos_sun(self):
        """
        Note: gLAB uses simplified expressions for determination of solar coordinates (based on Section 3.3.2 in
              :cite:`montenbruck2012`), whereas Where uses more sophisticated JPL ephemeris. The gLAB routine 'findsun'
              is implemented in Where and can be used for testing like

                glab_sun_earth = gnss.findsun(self.dset3.time)
                glab_pos_sun = glab_sun_earth - self.dset3.sat_posvel.itrs_pos
                glab_pos_sun_unit = mathp.unit_vector(glab_pos_sun)
        """
        itrs_pos_sun = self.dset3.sat_posvel._itrs_pos_sun
        expected_itrs_pos_sun = np.array([[-0.959930, 0.044658, -0.276657]])
        np.testing.assert_allclose(itrs_pos_sun, expected_itrs_pos_sun, rtol=0, atol=1e-2)
        # print('OUTPUT:\n itrs_pos_sun = {} m'.format(itrs_pos_sun))

    def test_gcrs_vel(self):
        """
        Note:  No external values used
        """
        gcrs_vel = self.dset3.sat_posvel.gcrs_vel
        expected_gcrs_vel = np.array([[-847.73680227, 760.19567755, 1.36311117]])
        np.testing.assert_allclose(gcrs_vel, expected_gcrs_vel, rtol=0, atol=1e-8)


if __name__ == "__main__":
    unittest.main()
