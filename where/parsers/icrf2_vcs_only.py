"""A parser for reading VCS radio source coordinates from ICRF2 files

Description:
------------

Reads radio source coordinates from ICRF2 :cite:`icrf2` files. ICRF2 is split into two sets of sources, namely VLBA
Calibrator Survey Sources (VCS sources) and non VCS sources. These two sets are listed in two separate files. The
defining sources are all non-VCS sources. In addition, some non-VCS sources are classified as special handling sources
because they exhibit temporal variations in celestial coordinates.

"""

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.parsers._parser_line import LineParser


@plugins.register
class Icrf2VscOnlyParser(LineParser):
    """A parser for reading source coordinates from ICRF files
    """

    def setup_parser(self):
        return dict(
            delimiter=(4, 17, 10, 4, 3, 12, 5, 3, 11, 12, 10, 8, 9, 8, 8, 7, 7),
            dtype="U5, U18, U10, i8, i8, f8, f8, i8, f8, f8, f8, b, f8, f8, f8, i8, i8",
            names=[
                "icrf",
                "icrf_name",
                "iers_name",
                "ra_h",
                "ra_m",
                "ra_s",
                "dec_deg",
                "dec_m",
                "dec_s",
                "ra_err",
                "dec_err",
                "special",
                "mean_mjd",
                "first_mjd",
                "last_mjd",
                "no_sessions",
                "no_obs",
            ],
            usecols=(1, 2, 3, 4, 5, 6, 7, 8),
            skip_header=20,
            autostrip=True,
        )

    def structure_data(self):
        ra = Unit.hms_to_rad(self._array["ra_h"], self._array["ra_m"], self._array["ra_s"])
        dec = Unit.dms_to_rad(self._array["dec_deg"], self._array["dec_m"], self._array["dec_s"])

        src_type = dict(vcs=True, non_vcs=False, defining=False, special=False, undefined=False)
        self.data = {
            src["iers_name"]: dict(icrf_name=src["icrf_name"], ra=ra[i], dec=dec[i], **src_type)
            for i, src in enumerate(self._array)
        }
