""" Test :mod:`where.apriori.orbit.broadcast`.

-------


$Revision: 18265 $
$Date: 2019-09-20 12:52:36 +0200 (Fri, 20 Sep 2019) $
$LastChangedBy: kirann $

"""
# Standard library imports
from datetime import datetime
import unittest

# External library imports
import numpy as np
import pytest

# Where imports
from where import apriori
from where.data.time import Time

TEST = "test_2"


@pytest.mark.quick
class TestBroadcast(unittest.TestCase):
    def setUp(self):
        """

        The first test set up is based on the bc_velo.c program, which is published in :cite:`remondi2004` and
         following RINEX navigation file sample:

        /* Sample Broadcast Message in unit of radians, seconds, meters.
        20 01  7 23  2  0  0.0 -.857324339449D-04 -.272848410532D-11  .000000000000D+00
             .200000000000D+02  .886875000000D+02  .465376527657D-08  .105827953357D+01
             .457651913166D-05  .223578442819D-02  .177137553692D-05  .515379589081D+04
             .936000000000D+05  .651925802231D-07  .164046615454D+01 -.856816768646D-07
             .961685061380D+00  .344968750000D+03  .206374037770D+01 -.856928551657D-08
             .342514267094D-09  .000000000000D+00  .112400000000D+04  .000000000000D+00
             .200000000000D+01  .000000000000D+00 -.651925802231D-08  .276000000000D+03
             .865800000000D+05  .000000000000D+00  .000000000000D+00  .000000000000D+00
        */


        The second test set up compares results from Where against gLAB solution for satellite G20 and epoch
        2016-03-01 00:00:00.0.

        /* Sample Broadcast Message in unit of radians, seconds, meters for satellite G20 and
        /  epoch 2016-03-01 00:00:00.0
        20 16  3  1  0  0  0.0 0.396233052015D-03 0.261479726760D-11 0.000000000000D+00
            0.100000000000D+02-0.231562500000D+02 0.530236372187D-08 0.253477496869D+00
           -0.111199915409D-05 0.483385741245D-02 0.810064375401D-05 0.515369705963D+04
            0.172800000000D+06-0.141561031342D-06 0.304306271006D+00 0.372529029846D-08
            0.926615731710D+00 0.207250000000D+03 0.133849764271D+01-0.843427989304D-08
           -0.164292557730D-09 0.100000000000D+01 0.188600000000D+04 0.000000000000D+00
            0.200000000000D+01 0.000000000000D+00-0.838190317154D-08 0.100000000000D+02
            0.172770000000D+06 0.400000000000D+01 0.000000000000D+00 0.000000000000D+00

        The second test set up compares results from Where against gLAB solution for satellite E11 and epoch
        2016-03-01 00:00:00.0.

        /* Sample Broadcast Message in unit of radians, seconds, meters for satellite E11 and
        /  epoch 2016-03-01 00:00:00.0
        E11 2016 03 01 00 00 00 6.643886445090e-05 1.097077984014e-11 0.000000000000e+00
             3.200000000000e+01-3.243750000000e+01 3.015839907552e-09 2.397505637802e+00
            -1.462176442146e-06 3.306962316856e-04 8.240342140198e-06 5.440621692657e+03
             1.728000000000e+05 9.313225746155e-09-1.259905024101e+00 5.960464477539e-08
             9.679475503522e-01 1.663125000000e+02-6.590241211713e-01-5.572732126661e-09
             2.775115594704e-10 2.580000000000e+02 1.886000000000e+03                   
             3.120000000000e+00 0.000000000000e+00-2.328306436539e-08 0.000000000000e+00
             1.735000000000e+05 

        """

        # Get GNSS ephemeris data for testing
        if TEST == "test_1":
            file_key = "test_apriori_orbit_broadcast_1"
            year = 2001
            month = 7
            day = 23
            hour = 2
            minute = 0
            second = 0.0
            satellite = "G20"
            self.system = "G"  # GNSS identifier

            # Satellite transmission time
            self.t_sat = 86400.00

        elif TEST == "test_2":
            file_key = "test_apriori_orbit_broadcast_2"
            year = 2016
            month = 3
            day = 1
            hour = 0
            minute = 0
            second = 0.0
            satellite = "G20"
            self.system = "G"  # GNSS identifier

            # Satellite transmission time
            self.t_sat = 172799.92312317

        elif TEST == "test_3":
            file_key = "test_apriori_orbit_broadcast_3"
            year = 2016
            month = 3
            day = 1
            hour = 0
            minute = 0
            second = 0.0
            satellite = "E11"
            self.system = "E"  # GNSS identifier

            # Satellite transmission time
            self.t_sat = 173699.999

        rundate = datetime(year, month, day, hour, minute)
        time = Time(
            [
                (
                    "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
                    "".format(year=year, month=month, day=day, hour=hour, minute=minute, second=second)
                )
            ],
            fmt="isot",
            scale="gps",
        )

        self.brdc = apriori.get(
            "orbit",
            apriori_orbit="broadcast",
            rundate=rundate,
            time=time,
            satellite=tuple({satellite}),
            system=tuple({self.system}),
            station="test",
            file_key=file_key,
        )

        self.idx = 0  # Broadcast ephemeris index

    def test_get_corrected_broadcast_ephemeris(self):
        """
        The test is based on the bc_velo.c program, which is published in :cite:`remondi2004`.
        """
        brdc_dict = self.brdc._get_corrected_broadcast_ephemeris(self.t_sat, self.idx, self.system)

        if TEST == "test_1":
            expected_a = np.array([26561612.084130041])  # Semimajor axis
            expected_E = np.array([8.190663e-03])  # Eccentric anomaly
            expected_i = np.array([9.616826e-01])  # Inclination
            expected_lambda = np.array([-4.659860e00])  # Instantaneous Greenwich longitude of the ascending node
            expected_n = np.array([1.458482e-04])  # Corrected mean motion
            expected_r = np.array([26501967.581685])  # Orbit radius
            expected_tk = np.array([-7.200000e03])  # Eclapsed time referred to ephemeris reference epoch
            expected_u = np.array([2.071945e00])  # Argument of latitude
            expected_vega = np.array([8.208996e-03])  # True anomaly

        elif TEST == "test_2":
            expected_a = np.array([26560593.382438906])  # Semimajor axis (a^2)
            expected_E = np.array([0.2546841246])  # Eccentric anomaly (E)
            expected_i = np.array([0.926616])  # Inclination (i)
            expected_lambda = np.array(
                [-12.296463]
            )  # Instantaneous Greenwich longitude of the ascending node (taken from Where solution) (O)
            expected_n = np.array([0.000145857])  # Corrected mean motion (calculated based on gLAB values)
            expected_r = np.array([26436138.824729])  # Orbit radius (r)
            expected_tk = np.array([-0.076876828621607])  # Eclapsed time referred to ephemeris reference epoch (diff)
            expected_u = np.array([1.594403])  # Argument of latitude (u)
            expected_vega = np.array([0.2559048275])  # True anomaly (fk)

        elif TEST == "test_3":
            expected_a = np.array([29600364.402609915])  # Semimajor axis (a^2)
            expected_E = np.array([2.509278395268e00])  # Eccentric anomaly (E)
            expected_i = np.array([0.967948])  # Inclination (i)
            expected_lambda = np.array(
                [-1.392631397645e01]
            )  # Instantaneous Greenwich longitude of the ascending node (O)
            expected_n = np.array([1.239749284615e-04])  # Corrected mean motion (n)
            expected_r = np.array([29608136.838406])  # Orbit radius (r)
            expected_tk = np.array([899.999000000000024])  # Eclapsed time referred to ephemeris reference epoch (diff)
            expected_u = np.array([1.850446560922e00])  # Argument of latitude (u)
            expected_vega = np.array([2.509473815034e00])  # True anomaly (fk)

        with self.subTest(msg="a"):
            np.testing.assert_allclose(brdc_dict["a"], expected_a, rtol=0, atol=1e-6)

        with self.subTest(msg="E"):
            np.testing.assert_allclose(brdc_dict["E"], expected_E, rtol=0, atol=1e-6)

        with self.subTest(msg="i"):
            np.testing.assert_allclose(brdc_dict["i"], expected_i, rtol=0, atol=1e-6)

        with self.subTest(msg="lambda_"):
            np.testing.assert_allclose(brdc_dict["lambda_"], expected_lambda, rtol=0, atol=1e-6)

        with self.subTest(msg="n"):
            np.testing.assert_allclose(brdc_dict["n"], expected_n, rtol=0, atol=1e-6)

        with self.subTest(msg="r"):
            np.testing.assert_allclose(brdc_dict["r"], expected_r, rtol=0, atol=1e-6)

        with self.subTest(msg="tk"):
            np.testing.assert_allclose(brdc_dict["tk"], expected_tk, rtol=0, atol=1e-6)

        with self.subTest(msg="u"):
            np.testing.assert_allclose(brdc_dict["u"], expected_u, rtol=0, atol=1e-6)

        with self.subTest(msg="vega"):
            np.testing.assert_allclose(brdc_dict["vega"], expected_vega, rtol=0, atol=1e-6)

        # print("OUTPUT:\n a = {:f} [m]\n E = {:15.12e} [rad]\n i = {:15.12e} [rad]\n lambda_ = {:15.12e} [rad]\n"
        #      " n = {:15.12e} [rad/s]\n r = {:f} [m]\n tk = {:15.12e} [s]\n u = {:15.12e} [rad]\n vega = {:15.12e} [rad]\n"
        #      ''.format(brdc_dict['a'], brdc_dict['E'], brdc_dict['i'], brdc_dict['lambda_'], brdc_dict['n'],
        #                brdc_dict['r'], brdc_dict['tk'], brdc_dict['u'], brdc_dict['vega']))

    def test_get_satellite_position_velocity(self):
        """
        The test is based on the bc_velo.c program, which is published in :cite:`remondi2004`.
        """
        sat_pos, sat_vel = self.brdc._get_satellite_position_velocity(self.t_sat, self.idx, self.system)

        if TEST == "test_1":
            expected_sat_pos = np.array([-12611434.19782218677, -13413103.97797041245, 19062913.07357876940])
            expected_sat_vel = np.array([266.2803795674, -2424.7683468482, -1529.7620784616])

        elif TEST == "test_2":
            expected_sat_pos = np.array([-4833470.67769, 15130611.05029, 21132240.92266])
            expected_sat_vel = np.array([-2632.8593, -735.4519, -69.6660])

        elif TEST == "test_3":
            expected_sat_pos = np.array([14067596.917149, 11368467.367032, 23441468.250584])
            expected_sat_vel = np.array([0, 0, 0])  # TODO: Not known

        with self.subTest(msg="sat_pos"):
            np.testing.assert_allclose(sat_pos, expected_sat_pos, rtol=0, atol=1e-5)

        with self.subTest(msg="sat_vel"):
            np.testing.assert_allclose(sat_vel, expected_sat_vel, rtol=0, atol=1e-6)

        # print("OUTPUT:\n sat_pos = {:20.8f} {:20.8f} {:20.8f} [m]\n sat_vel = {:15.8f} {:15.8f} {:15.8f} [m/s]\n"
        #      ''.format(sat_pos[0], sat_pos[1], sat_pos[2], sat_vel[0], sat_vel[1], sat_vel[2]))

    @pytest.mark.sisre
    @pytest.mark.gnss
    # @pytest.mark.xfail(reason="Tada")
    def test_get_satellite_clock_correction(self):

        sat_clk_corr = self.brdc.satellite_clock_correction()

        if TEST == "test_1":
            expected_sat_clk_corr = -25696.047654177062
        elif TEST == "test_2":
            expected_sat_clk_corr = 118787.68054
        elif TEST == "test_3":
            expected_sat_clk_corr = 19920.830541813138552

        np.testing.assert_allclose(sat_clk_corr, expected_sat_clk_corr, rtol=0, atol=1e-5)

    def test_get_relativistic_clock_correction(self):
        rel_clk_corr = self.brdc._get_relativistic_clock_correction(self.t_sat, self.idx, self.system)

        if TEST == "test_1":
            expected_rel_clk_corr = -0.012570413426601913
        elif TEST == "test_2":
            expected_rel_clk_corr = -0.83775
        elif TEST == "test_3":
            expected_rel_clk_corr = 0  # TODO: Not known

        np.testing.assert_allclose(rel_clk_corr, expected_rel_clk_corr, rtol=0, atol=1e-6)


if __name__ == "__main__":
    unittest.main()
