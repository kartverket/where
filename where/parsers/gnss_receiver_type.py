"""A parser for reading GNSS receiver type file

Example:
--------

    from where import parsers
    parser = parsers.parse('gnss_receiver_type')

Description:
------------

Reads data from GNSS receiver type file.


"""

# Standard library imports
import itertools
import re

# Where imports
from where.parsers import parser
from where.lib import plugins


@plugins.register
class GnssReceiverTypeParser(parser.Parser):
    """
    """

    def __init__(self):
        super().__init__()
        self.file_key = "gnss_receiver_types"

    #
    # PARSERS
    #
    def setup_parsers(self):
        """Parsers defined for reading GNSS receiver type file line by line.
        """
        file_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: not re.match("^#|^\s*$", line),  # skip comments starting with '#' and empty lines
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9
                # ALTUS APS-3          P  Y     # GPS L1/L2+L2C, GLO L1/L2, SBAS integrated rcvr/antenna
                # ALTUS APS-3L         P  Y     # GPS L1/L2+L2C, GLO L1/L2, SBAS, L-Band rcvr/antenna
                True: {
                    "parser": self.parse_string, "fields": {"name": (0, 20), "type": (20, 22), "igs_format": (22, 25)}
                }
            },
        )

        return itertools.repeat(file_parser)

    def parse_string(self, line, _):
        """Parse string entries of GNSS receiver type file to instance variable 'data'.
        """
        self.data.update({line["name"]: {"type": line["type"], "igs_format": line["igs_format"]}})

    #
    # WRITE DATA
    #
    def write_data(self, dout):
        """Write data based on GNSS receiver type file

        Args:
            dout (dict): Dictionary with GNSS receiver name as key and following entries as values:

        =============  ================================================================================================
         Value          Description
        =============  ================================================================================================
         igs_format     GNSS receiver name follows IGS convention indicated with Y(es) or N(ot). The IGS receiver name
                        convention is defined in IGS "rcvr_ant.tab" file.
         type           GNSS receiver type indicated by 'X', 'C' or 'P', which stands for:
                           'X' : receiver is cross-correlating and requires correction of P2' and C1 Rogue SNR,
                                 Trimble 4000, etc.
                           'C' : receiver is non-cross-correlating but reports C1 instead of P1 Trimble 4700, 5700,
                                 Leica RS500, CRS1000, SR9600, etc. unless AS is off
                           'P' : receiver is non-cross-correlating and reports true P1, P2
        =============  ================================================================================================
        """
        dout.update(self.data)
