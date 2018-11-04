"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads epoch of discontinuity for the position and velocity model in ITRF.

The file is using the SINEX format, but the SOLUTION/DISCONTINUITY block is not defined in the official format
description.

"""

# Standard library imports
from datetime import datetime

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_sinex import SinexParser, SinexBlock, SinexField


@plugins.register
class ItrfSnxSolnParser(SinexParser):
    """A parser for reading data from ITRF files in SNX format
    """

    def __init__(self, file_path, encoding=None, logger=None, header=None):
        """Set up the basic information needed by the parser

        Turn off parsing of header by default, as the soln-Sinex files in general are missing headers.

        Args:
            file_path (String/Path):    Path to file that will be read.
            encoding (String):          Encoding of file that will be read.
            logger (Function):          Ignored for where.parsers, used for consistency with Midgard.
            header (Boolean):           Whether to parse the header.
        """
        super().__init__(file_path, encoding)
        self._header = False if header is None else header

    def setup_parser(self):
        return {self.solution_discontinuity}

    @property
    def solution_discontinuity(self):
        """Custom made block for ITRF to mark the epoch of discontinuities in the position and velocity of stations

        Example:
             1515  A    1 R 00:000:00000 92:180:43054 P - EQ M7.3 - Southern California
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
        """
        return SinexBlock(
            marker="SOLUTION/DISCONTINUITY",
            fields=(
                SinexField("site_code", 1, "U4"),
                SinexField("point_code", 6, "U2"),
                SinexField("soln", 9, "U4"),
                SinexField("obs_code", 14, "U1"),
                SinexField("start_epoch", 16, "O", "epoch"),
                SinexField("end_epoch", 29, "O", "epoch"),
                SinexField("pos_or_vel", 42, "U1"),
                SinexField("_dash", 44, None),  # Ignored
                SinexField("comment", 46, "U34"),
            ),
            parser=self.parse_solution_discontinuity,
        )

    def parse_solution_discontinuity(self, data):
        """Parser for SOLUTION/DISCONTINUITY data

        Converts the input data to a dictionary with items for each site, each containing start and end epochs for each
        solution. Only uses the position data, not the velocity data.

        Args:
            data (numpy.array):  Input data, raw data for SOLUTION/DISCONTINUITY block.
        """
        pos_data = data[data["pos_or_vel"] == "P"]
        pos_data["start_epoch"][np.equal(pos_data["start_epoch"], None)] = datetime.min
        pos_data["end_epoch"][np.equal(pos_data["end_epoch"], None)] = datetime.max

        for d, soln in zip(pos_data, pos_data["soln"].astype("i8")):
            site_key = d["site_code"]
            self.data.setdefault(site_key, dict())
            self.data[site_key].setdefault(soln, dict(start=d["start_epoch"], end=d["end_epoch"]))
