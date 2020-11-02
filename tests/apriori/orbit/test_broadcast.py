""" Test :mod:`where.apriori.orbit.broadcast`.

-------


$Revision: 19696 $
$Date: 2020-10-26 22:48:44 +0100 (ma., 26 okt. 2020) $
$LastChangedBy: dahmic $

"""
# Standard library imports
from datetime import datetime
import pathlib
import unittest

# External library imports
import numpy as np
import pytest

# Midgard imports
from midgard.math.constant import constant

# Where imports
from where import apriori
from where.data import dataset3 as dataset
from where.lib import config
from where.lib import log

TEST = "test_3"


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

        The third test set up compares results from Where against CNES solution for satellite E01 and epoch
        2016-03-01 01:30:00.0.

        /* Sample Broadcast Message in unit of radians, seconds, meters for satellite E01 and
        /  epoch 2019-07-01 00:00:00.0
        E01 2019 07 01 00 00 00-6.374700460583D-04-8.085976332950D-12 0.000000000000D+00
             1.600000000000D+01 2.106250000000D+02 2.413671967697D-09 5.641607228729D-01
             9.929761290550D-06 1.870252890512D-04 7.383525371552D-06 5.440612319946D+03
             8.640000000000D+04-3.725290298462D-09 2.424721177655D-01 1.657754182816D-07
             9.878562635157D-01 1.958125000000D+02 3.073143357419D+00-5.291648989849D-09
             2.003654888840D-10 2.580000000000D+02 2.060000000000D+03 0.000000000000D+00
             3.120000000000D+00 0.000000000000D+00-1.862645149231D-09 0.000000000000D+00
             8.714000000000D+04 0.000000000000D+00 0.000000000000D+00 0.000000000000D+00 

        """
        # Initialize logging
        log.init(log_level="debug")

        # Get GNSS ephemeris data for testing
        if TEST == "test_1":
            file_name = "test2040.01n"
            year = 2001
            month = 7
            day = 23
            hour = 2
            minute = 0
            second = 0
            satellite = "G20"
            self.system = "G"  # GNSS identifier

            # Satellite transmission time
            self.t_sat_gpsweek = 1124.0
            self.t_sat_gpssec = 86400.00

        elif TEST == "test_2":
            file_name = "test0610.16n"
            year = 2016
            month = 3
            day = 1
            hour = 0
            minute = 0
            second = 0
            satellite = "G20"
            self.system = "G"  # GNSS identifier

            # Satellite transmission time
            self.t_sat_gpsweek = 1886.0
            self.t_sat_gpssec = 172799.92312317

        elif TEST == "test_3":
            file_name = "TEST00CNS_R_20191820000_01D_EN.rnx"
            year = 2019
            month = 7
            day = 1
            hour = 1
            minute = 30
            second = 0
            satellite = "E01"
            self.system = "E"  # GNSS identifier

            # Satellite transmission time
            self.t_sat_gpsweek = 2060.0
            self.t_sat_gpssec = 91800.0

        rundate = datetime(year, month, day, hour, minute, second)

        # Initialize configuration
        config.init(rundate=rundate, pipeline="gnss")

        # Generate observation datast
        self.dset = dataset.Dataset(num_obs=1, rundate=rundate)
        self.dset.add_time(name="time", val=rundate, scale="gps", fmt="datetime")
        #self.dset.add_time(
        #    name="time", val=Time(val=[self.t_sat_gpsweek], val2=[self.t_sat_gpssec], fmt="gps_ws", scale="gps")
        #)
        self.dset.add_text(name="satellite", val=[satellite])

        # Get broadcast ephemeris
        self.brdc = apriori.get(
            "orbit",
            rundate=rundate,
            system=tuple({self.system}),
            station="test",
            apriori_orbit="broadcast",
            file_path=pathlib.Path(__file__).parent / "files" / file_name,
        )

        self.idx = 0  # Broadcast ephemeris index

    def test_get_corrected_broadcast_ephemeris(self):
        """
        The test is based on the bc_velo.c program, which is published in :cite:`remondi2004`.
        """
        brdc_dict = self.brdc._get_corrected_broadcast_ephemeris(
            self.t_sat_gpsweek, self.t_sat_gpssec, self.idx, self.system
        )

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
            expected_a = np.array([29600262.415948])  # Semimajor axis (a^2)
            expected_E = np.array([1.233802050337e00])  # Eccentric anomaly (E)
            expected_i = np.array([9.878574681865e-01])  # Inclination (i)
            expected_lambda = np.array(
                [-6.451718161810]
            )  # Instantaneous Greenwich longitude of the ascending node (O)
            expected_n = np.array([1.239725533334e-04])  # Corrected mean motion (n)
            expected_r = np.array([2.9598449621537e07])  # Orbit radius (r)
            expected_tk = np.array([5400.0])  # Eclapsed time referred to ephemeris reference epoch (diff)
            expected_u = np.array([4.30712042665])  # Argument of latitude (u)
            expected_vega = np.array([1.233978561432])  # True anomaly (fk)

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

        log.debug(f"\nOUTPUT test_get_corrected_broadcast_ephemeris():\n"
                  f"a =       {brdc_dict['a']:f} [m]\n"
                  f"E =       {brdc_dict['E']:15.12e} [rad]\n"
                  f"i =       {brdc_dict['i']:15.12e} [rad]\n"
                  f"lambda_ = {brdc_dict['lambda_']:15.12e} [rad]\n"
                  f"n =       {brdc_dict['n']:15.12e} [rad/s]\n" 
                  f"r =       {brdc_dict['r']:f} [m]\n"
                  f"tk =      {brdc_dict['tk']:f} [s]\n"
                  f"u =       {brdc_dict['u']:15.12e} [rad]\n"
                  f"vega =    {brdc_dict['vega']:15.12e} [rad]\n"
        )



    def test_get_satellite_position_velocity(self):
        """
        The test is based on the bc_velo.c program, which is published in :cite:`remondi2004`.
        """
        sat_pos, sat_vel = self.brdc._get_satellite_position_velocity(
            self.t_sat_gpsweek, self.t_sat_gpssec, self.idx, self.system
        )

        if TEST == "test_1":
            expected_sat_pos = np.array([-12611434.19782218677, -13413103.97797041245, 19062913.07357876940])
            expected_sat_vel = np.array([266.2803795674, -2424.7683468482, -1529.7620784616])

        elif TEST == "test_2":
            expected_sat_pos = np.array([-4833470.67769, 15130611.05029, 21132240.92266])
            expected_sat_vel = np.array([-2632.8593, -735.4519, -69.6660])

        elif TEST == "test_3":
            expected_sat_pos = np.array([-14015916.0717073, -12803962.9429849, -22708607.3906834])
            expected_sat_vel = np.array([0, 0, 0])  # TODO: Not known

        with self.subTest(msg="sat_pos"):
            np.testing.assert_allclose(sat_pos, expected_sat_pos, rtol=0, atol=1e-5)

        if not TEST == "test_3":
            with self.subTest(msg="sat_vel"):
                np.testing.assert_allclose(sat_vel, expected_sat_vel, rtol=0, atol=1e-6)

        log.debug(f"OUTPUT test_get_satellite_position_velocity():\n"
                  f"sat_pos = {sat_pos[0]:20.8f} {sat_pos[1]:20.8f} {sat_pos[2]:20.8f} [m]\n"
                  f"sat_vel = {sat_vel[0]:20.8f} {sat_vel[1]:20.8f} {sat_vel[2]:20.8f} [m/s]\n"
        )

    @pytest.mark.sisre
    @pytest.mark.gnss
    # @pytest.mark.xfail(reason="Tada")
    def test_get_satellite_clock_correction(self):

        sat_clk_corr = self.brdc.satellite_clock_correction(self.dset)

        if TEST == "test_1":
            expected_sat_clk_corr = -25696.047654177062
        elif TEST == "test_2":
            expected_sat_clk_corr = 118787.68054
        elif TEST == "test_3":
            expected_sat_clk_corr = -6.375137103305e-04 * constant.c # -191121.8022286806 -> conversion from s to m

        np.testing.assert_allclose(sat_clk_corr, expected_sat_clk_corr, rtol=0, atol=1e-5)
        
        log.debug(f"OUTPUT test_get_satellite_clock_correction():\n"
                  f"sat_clk_corr = {sat_clk_corr[0]:12.10f} [m]\n"
        )

    def test_get_relativistic_clock_correction(self):
        rel_clk_corr = self.brdc._get_relativistic_clock_correction(
            self.t_sat_gpsweek, self.t_sat_gpssec, self.idx, self.system
        )

        if TEST == "test_1":
            expected_rel_clk_corr = -0.012570413426601913
        elif TEST == "test_2":
            expected_rel_clk_corr = -0.83775
        elif TEST == "test_3":
            expected_rel_clk_corr = -4.266422241696e-10 * constant.c # = -0.1279041211 -> conversion from s to m

        np.testing.assert_allclose(rel_clk_corr, expected_rel_clk_corr, rtol=0, atol=1e-6)

        log.debug(f"OUTPUT test_get_relativistic_clock_correction():\n"
                  f"rel_clk_corr = {rel_clk_corr:12.10f} [m]\n"
        )


if __name__ == "__main__":
    unittest.main()
