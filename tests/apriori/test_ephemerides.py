""" Test :mod:`where.apriori.ephemerides`.


"""

# Standard library imports
import unittest

# External library imports
import numpy as np

# Where imports
from where import apriori
from where.data.time import Time


class TestEphemerides(unittest.TestCase):
    def setUp(self):

        # TODO: Configuration has to be defined? How? where.set_config(2016, 3, 1, 'gps')?
        time = Time(2457448.5, fmt="jd", scale="tdb")  # Julian Day in TDB time scale
        self.eph = apriori.get("ephemerides", time=time, ephemerides="de430")

    def test_ephemerides(self):
        """Test the use of JPL DE ephemeris with the `jplephem` package

        The vector from the Earth mass center to the Sun mass center in kilometers is determined in the BCRS for the
        test. This is done by carrying out following steps:
            * vector from Earth center of mass (NAIF ID: 399) to the Earth-Moon barycenter (NAIF ID: 3)
            * vector from Earth-Moon barycenter (NAIF ID: 3) to Solar System barycenter (NAIF ID: 0)
            * vector from Solar System barycenter (NAIF ID: 0) to Sun center of mass (NAIF ID: 10)

        The result is compared against results from HORIZONS Web-Interface
        (http://ssd.jpl.nasa.gov/horizons.cgi#results) for DE431 ephemeris with following settings:
            Ephemeris Type [change] :  	    VECTORS
            Target Body [change] :  	    Sun [Sol] [10]
            Coordinate Origin [change] :    Geocentric [500]
            Time Span [change] :            Start=2016-03-01, Stop=2016-03-02, Step=1 m
            Table Settings [change] :       output units=KM-S; reference plane=FRAME
            Display/Output [change] :  	    default (formatted HTML)

        Result:

            2457448.500000000 = A.D. 2016-Mar-01 00:00:00.0000 TDB                    # Julian Day in TDB time scale
               1.398250202951112E+08 -4.514599833424634E+07 -1.957150223235953E+07    # X, Y, Z in [km]
               1.038335388022904E+01  2.587894057258682E+01  1.121912358196623E+01    # VX, VY, VZ in [km/s]

        """

        # Determine Sun position vector in BCRS pointing from Earth center of mass to Sun center of mass
        sun_pos_bcrs = self.eph.bcrs("earth", "sun")
        expected_sun_pos_bcrs = np.array([1.398250202951112e11, -4.514599833424634e10, -1.957150223235953e10])

        # NOTE: Where can use other JPL ephemeris (e.g. DE430) instead of test ephemeris DE431, so some differences
        #       should be expected.
        np.testing.assert_allclose(sun_pos_bcrs, expected_sun_pos_bcrs, rtol=0, atol=1e-4)

        # print('OUTPUT:\n sun_pos_bcrs = {:f} [km] {:f} [km] {:f} [km]\n'
        #       ''.format(sun_pos_bcrs[0], sun_pos_bcrs[1], sun_pos_bcrs[2]))


if __name__ == "__main__":
    unittest.main()
