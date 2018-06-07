"""A parser for reading GNSS clock data from RINEX clock format

Example:
--------

    from where import parsers
    parser = parsers.parse('rinex_clk', rundate=rundate)

Description:
------------

The IGS RINEX clock file includes clock bias for IGS stations and GPS satellites aligned to GPS time. The format
description can be found in :cite:`ray2010`.

This routine reads so far only GPS satellite clock bias values given for every 5 minutes. The satellite clock bias is
given in seconds, but is converted to meters.



"""

# Standard library imports
from datetime import datetime
import itertools

# External library imports
import numpy as np

# Where imports
from where.parsers import parser
from where.lib import constant
from where.lib import plugins


@plugins.register
class RinexClkParser(parser.Parser):
    """A parser for reading GNSS satellite clock data from RINEX clock files
    """

    def __init__(self, rundate):
        """
        Args:
            rundate:       The model run date.
        """
        super().__init__(rundate)
        self.file_key = "gnss_rinex_clk"

    #
    # PARSER
    #
    def setup_parsers(self):
        """Define parsers for reading RINEX clock file
        """
        data_parser = parser.define_parser(
            end_marker=lambda _l, _ln, n: n.startswith("*"),
            label=lambda line, _ln: line[0:3],
            parser_def={
                "AS ": {
                    "parser": self.parse_sat_clk,
                    "fields": [
                        None,
                        "satellite",
                        "year",
                        "month",
                        "day",
                        "hour",
                        "minute",
                        "seconds",
                        "num_val",
                        "sat_clock_bias",
                        "sat_clock_bias_sigma",
                    ],
                }
            },
        )

        return itertools.repeat(data_parser)

    def parse_sat_clk(self, line, _):
        """Parse satellite clock entries marked with 'AS' in RINEX clock file

        Args:
            line (dict):  Dict containing the fields of a line.
        """
        sat = line["satellite"].upper()
        time = datetime(
            int(line["year"]),
            int(line["month"]),
            int(line["day"]),
            int(line["hour"]),
            int(line["minute"]),
            round(float(line["seconds"])),
        )

        self.data.setdefault("time", list()).append(time)
        self.data.setdefault("satellite", list()).append(sat)
        self.data.setdefault("sat_clock_bias", list()).append(
            float(line["sat_clock_bias"]) * constant.c
        )  # Convert [s]
        # to [m]

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Write data based on GNSS satellite clock data from RINEX clock files

        Args:
            dset (Dataset): Dataset with precise satellite clock values stored as follows

        ====================  ===============  =======  ========================================================
         Field                 Type             Unit     Description
        ====================  ===============  =======  ========================================================
         sat_clock_bias        numpy.ndarray     m       Satellite clock offset from GPS time
         satellite             numpy.ndarray             Satellite PRN number
         time                  TimeTable                 Observation epoch in GPS time
        ====================  ===============  =======  ========================================================

        """

        dset.num_obs = len(self.data["time"])

        for k, v in self.data.items():
            if k == "time":
                dset.add_time(k, val=v, scale="gps")
            elif k == "satellite":
                dset.add_text(k, val=v)
            else:
                if isinstance(v, list):
                    dset.add_float(k, val=np.array(v))
