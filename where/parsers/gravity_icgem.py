"""A parser for reading coefficients of the spherical harmonic functions which determine the Earth gravitational field.

Description:
------------

asdf


References:
-----------

http://icgem.gfz-potsdam.de/ICGEM/documents/ICGEM-Format-2011.pdf
http://icgem.gfz-potsdam.de/ICGEM/


http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/

Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
                IERS Technical Note No. 36, BKG (2010)

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log
from where.parsers import parser


@plugins.register
class GravityIcgemParser(parser.ParserDict):
    """A parser for reading gravity coefficients from files
    """

    def __init__(self, gravity_field, truncation_level):
        super().__init__()

        self.vars["gravity_field"] = gravity_field
        self.truncation_level = truncation_level
        self.data["C"] = np.zeros((truncation_level + 1, truncation_level + 1))
        self.data["S"] = np.zeros((truncation_level + 1, truncation_level + 1))
        self.data["C"][0, 0] = 1  # C[0, 0] is always 1, but not always given on file

        # Keep track of which coefficients that have been read
        self.coeff_read = np.triu(np.ones((truncation_level + 1, truncation_level + 1)), 1)
        self.coeff_read[0:2, :] = True  # Coeffs of degree 0 and 1 not on file

    #
    # PARSER for reading geopotential coefficients
    #
    def setup_parsers(self):
        header_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line.startswith("end_of_head"),
            label=lambda line, _ln: (line + " .").split()[0],
            parser_def={
                "earth_gravity_constant": {"parser": self.parse_constant, "fields": [None, "GM"]},
                "radius": {"parser": self.parse_constant, "fields": [None, "a"]},
            },
        )

        coefficient_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: np.all(self.coeff_read),
            label=lambda line, _ln: (line + " .").split()[0],
            parser_def={
                "gfc": {"parser": self.parse_coeff, "fields": [None, "degree", "order", "coeff_C", "coeff_S"]}
            },
        )

        return itertools.chain([header_parser, coefficient_parser])

    def parse_coeff(self, line, _):
        """Parse one line of gravitational data
        """
        n = int(line["degree"])
        m = int(line["order"])

        if n > self.truncation_level or m > self.truncation_level:
            return

        # TODO: Conversion from d to e only needed for C_00 in egm2008.gfc?
        self.data["C"][n, m] = float(line["coeff_C"].replace("d", "e"))
        self.data["S"][n, m] = float(line["coeff_S"].replace("d", "e"))
        self.coeff_read[n, m] = True

    def parse_constant(self, line, _):
        log.warn("TODO: parse_constant {}", line)
