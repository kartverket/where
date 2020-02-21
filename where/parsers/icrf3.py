"""A parser for reading radio source coordinates from ICRF3 files

Description:
------------

"""

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.parsers._parser_line import LineParser


@plugins.register
class Icrf3Parser(LineParser):
    """A parser for reading source coordinates from ICRF3 files
    """

    def setup_parser(self):
        return dict(
            delimiter=(4, 17, 12, 3, 6, 3, 12, 7, 3, 11, 15, 14, 8, 10, 9, 9, 6, 7, 7),
            dtype="U5, U17, U12, b, i8, i8, f8, f8, i8, f8, f8, f8, b, f8, f8, f8, i8, i8, i8",
            names=[
                "icrf",
                "icrf_name",
                "iers_name",
                "defining",
                "ra_h",
                "ra_m",
                "ra_s",
                "dec_deg",
                "dec_m",
                "dec_s",
                "ra_err",
                "dec_err",
                "corr_ra_dec",
                "mean_mjd",
                "first_mjd",
                "last_mjd",
                "no_sessions",
                "no_obs",
                "unknown",
            ],
            converters={
                3: lambda s: True if s.decode("utf-8").strip() == "D" else False,
                # 12: lambda s: True if s.decode("utf-8").strip() == "0.000" else False,
            },
            usecols=(1, 2, 3, 4, 5, 6, 7, 8, 9),
            skip_header=22,
            autostrip=True,
        )

    def structure_data(self):
        ra = Unit.hms_to_rad(self._array["ra_h"], self._array["ra_m"], self._array["ra_s"])
        dec = Unit.dms_to_rad(self._array["dec_deg"], self._array["dec_m"], self._array["dec_s"])

        src_type = dict(vcs=False, non_vcs=False, undefined=False, special=False)
        self.data = {
            src["iers_name"]: dict(
                icrf_name=src["icrf_name"], defining=src["defining"], ra=ra[i], dec=dec[i], **src_type
            )
            for i, src in enumerate(self._array)
        }
