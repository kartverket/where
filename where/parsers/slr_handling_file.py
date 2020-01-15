"""A parser for reading SLR handling file

Description:
------------

A parser for reading the ILRS Data handling file from DGFI. This file summarizes information on range-, time- or pressure biases and data to be deleted. It is adopted by the ILRS Analysis Working Group. 

References:
-----------

http://ilrs.dgfi.tum.de/fileadmin/data_handling/ILRS_Data_Handling_File.snx

"""
# Standard library imports
from datetime import datetime

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_sinex import SinexParser, SinexBlock, SinexField
from midgard.dev import log


@plugins.register
class SlrHandlingFileParser(SinexParser):
    """A parser for reading SLR handling file
    """

    def __init__(self, file_path, encoding=None, header=None):
        """Set up the basic information needed by the parser

        Turn off parsing of header by default, as the header contains a mistake

        Args:
            file_path (String/Path):    Path to file that will be read.
            encoding (String):          Encoding of file that will be read.
            header (Boolean):           Whether to parse the header.
        """
        super().__init__(file_path, encoding)
        self._header = False if header is None else header

    def setup_parser(self):
        return (self.solution_data_handling,)

    @property
    def solution_data_handling(self):
        """

        Example:
            *CODE PT_ UNIT T _DATA_START_ __DATA_END__ M __E-VALUE___ STD_DEV ___COMMENTS______
             7080 --- mm   A 16:106:00000 00:000:00000 E                      small remaining range bias
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
        """
        return SinexBlock(
            marker="SOLUTION/DATA_HANDLING",
            fields=(
                SinexField("site_code", 1, "U4"),
                SinexField("point_code", 6, "U3"),  # TODO "---" means all satellites, not sure what L55 etc means.
                SinexField("unit", 10, "U4"),
                SinexField(
                    "obs_code", 15, "U1"
                ),  # TODO: Indicates if information involves first or second wavelength?
                SinexField("start_time", 17, "O", "epoch"),
                SinexField("end_time", 30, "O", "epoch"),
                SinexField("handling_code", 43, "U1"),  # Indicated what is to be done with the data.
                SinexField("e_value", 45, "U12"),
                SinexField("std_dev", 58, "U7"),
                SinexField("comments1", 66, "U9"),
                SinexField("comments2", 76, "U7"),
            ),
            parser=self.parse_data_handling,
        )

    def parse_data_handling(self, data):
        for d in data:
            start_time = datetime.min if d["start_time"] is None else d["start_time"]
            end_time = datetime.max if d["end_time"] is None else d["end_time"]
            interval = (start_time, end_time)
            info = {"unit": d["unit"]}

            if d["e_value"]:
                try:
                    e_value = float(d["e_value"])
                except ValueError:
                    log.fatal("ILRS Data handling: Not able to convert value to float")
                info.update({"e_value": e_value})
            if d["std_dev"]:
                try:
                    std_dev = float(d["std_dev"])
                except ValueError:
                    log.fatal("ILRS Data handling: Not able to convert value to float")
                info.update({"std_dev": std_dev})

            # Unfortunately we have to deal with two different line formats.
            # Split the comments field in the second line format:
            # *CODE PT_ UNIT T _DATA_START_ __DATA_END__ M __E-VALUE___ STD_DEV ___COMMENTS______
            # *CODE PT_ UNIT T _DATA_START_ __DATA_END__ M __E-VALUE___ STD_DEV _E-RATE__ _CMNTS_
            try:
                info.update({"comments": d["comments2"], "e_rate": float(d["comments1"])})
            except ValueError:
                info.update({"comments": d["comments1"] + d["comments2"]})

            self.data.setdefault(d["site_code"], {}).setdefault(d["handling_code"], []).append((interval, info))
