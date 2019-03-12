"""A module to provide information from radio sources defined in ICRF3

Description:
------------

Reads source positions ICRF3 files. Source positions are reported at epoch 2015.0 in ICRF3 and contains a model
for account for Galactocentric acceleration.

Galactic center position from wikipedia:
In 1958 the International Astronomical Union (IAU) decided to adopt the position of Sagittarius A as the true zero
co-ordinate point for the system of galactic latitude and longitude.[6] In the equatorial coordinate system the
location is: RA  17h 45m 40.04s, Dec −29° 00′ 28.1″ (J2000 epoch).

References:
-----------




"""

# External library imports
import numpy as np

# Midgard imports
from midgard.math.unit import Unit

# Where imports
from where.apriori import crf
from where.lib import plugins
from where import parsers

# TODO: Use dataset3 Direction and add derivatives as properties
def radec2unitvector(ra, dec):
    """Convert right ascension and declination to a unit vector
    """
    return np.array([np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)])


def dsrc_dra(ra, dec):
    """Derivative of radio source unit vector with regards to right ascension
    """
    cos_dec = np.cos(dec)
    cos_ra = np.cos(ra)
    # sin_dec = np.sin(dec)
    sin_ra = np.sin(ra)
    return np.array([-cos_dec * sin_ra, cos_dec * cos_ra, 0])


def dsrc_ddec(ra, dec):
    cos_dec = np.cos(dec)
    cos_ra = np.cos(ra)
    sin_dec = np.sin(dec)
    sin_ra = np.sin(ra)
    return np.array([-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec])


@plugins.register
class Icrf3(crf.CrfFactory):
    """A class to provide information from radio sources defined in ICRF2
    """

    def __init__(self, time, catalog=None):
        super().__init__(time)
        self.catalog = catalog

    def _read_data(self):
        """Read data needed by this Celestial Reference Frame for calculating positions of sources

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        data = parsers.parse_key(file_key="icrf3", file_vars=dict(catalog=self.catalog)).as_dict()

        return data

    def _calculate_pos_crs(self, source):
        """Calculate position for a source

        Args:
            source (String):    Key saying which source to calculate position for.

        Returns:
            Array:  Positions, one 2-vector
        """

        ga = 0.0058 * Unit.mas2rad  # Aberration constant (mas/yr)
        mjd_2015 = 57023.0  # # Reference epoch of aberration model

        # Galactic center
        gc_ra = Unit.hms_to_rad(17, 45, 40.04)
        gc_dec = Unit.dms_to_rad(-29, 0, 28.1)
        gc = radec2unitvector(gc_ra, gc_dec)

        source_info = self.data[source]

        # Compute correction
        dra = ga * gc @ dsrc_dra(source_info["ra"], source_info["dec"])
        ddec = ga * gc @ dsrc_ddec(source_info["ra"], source_info["dec"])
        dt = (self.time.mean.mjd - mjd_2015) * Unit.day2julian_year

        return np.array([source_info["ra"] + dra * dt, source_info["dec"] + ddec * dt])
