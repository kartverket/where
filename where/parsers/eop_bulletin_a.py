"""A parser for reading data from EOP files

Description:
------------

Reads data from EOP files.

"""

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where.parsers._parser_line import LineParser


@plugins.register
class EopBulletinAParser(LineParser):
    """A parser for reading data from EOP files
    """

    def setup_parser(self):

        return dict(
            delimiter=(
                2,
                2,
                2,
                1,
                8,
                1,
                1,
                1,
                9,
                9,
                1,
                9,
                9,
                2,
                1,
                10,
                10,
                1,
                7,
                7,
                2,
                1,
                1,
                9,
                9,
                1,
                9,
                9,
                10,
                10,
                11,
                10,
                10,
            ),
            dtype="i2, i2, i2, u1, f8, u1, u1, u1, f8, f8, u1, f8, f8, u2, u1, f8,f8, u1, f8,f8, u2,u1,u1,f8,f8,u1,f8,f8,f8,f8,f8,f8,f8",
            names=[
                "year",
                "month",
                "day",
                "blank",
                "mjd",
                "blank",
                "pm_flag",
                "blank",
                "x",
                "x_ferr",
                "blank",
                "y",
                "y_ferr",
                "blank",
                "ut1_utc_flag",
                "ut1_utc",
                "ut1_utc_ferr",
                "blank",
                "lod",
                "lod_ferr",
                "blank",
                "nut_flag",
                "blank",
                "dx",
                "dx_ferr",
                "blank",
                "dy",
                "dy_ferr",
                "b_x",
                "b_y",
                "b_ut1_utc",
                "b_dx",
                "b_dy",
            ],
            autostrip=True,
        )

    def structure_data(self):
        self._array["lod"] *= Unit.ms2s
        self._array["dx"] *= Unit.milliarcsec2arcsec
        self._array["dy"] *= Unit.milliarcsec2arcsec
        self.data = {
            item["mjd"]: dict(
                x=item["x"], y=item["y"], ut1_utc=item["ut1_utc"], lod=item["lod"], dx=item["dx"], dy=item["dy"]
            )
            for item in self._array
        }
