"""A parser for reading coefficients of the spherical harmonic functions which determine the Earth gravitational field.

Description:
------------

TODO: Should we use constants in header of file instead of where constants?


References:
-----------

http://icgem.gfz-potsdam.de/ICGEM-Format-2011.pdf
http://icgem.gfz-potsdam.de/tom_longtime

http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/

Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
                IERS Technical Note No. 36, BKG (2010)

"""

# Standard library imports
import itertools
import pathlib
from typing import Callable, NamedTuple, Optional, Union

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser


class HeaderField(NamedTuple):
    name: str
    converter: Callable


HEADER_FIELDS = (
    HeaderField("earth_gravity_constant", float),
    HeaderField("radius", float),
    HeaderField("max_degree", int),
)


@plugins.register
class GravityIcgemParser(Parser):
    """A parser for reading gravity coefficients from files
    """

    def __init__(
        self, file_path: Union[str, pathlib.Path], num_degrees: Optional[int] = None, encoding: Optional[str] = None
    ) -> None:
        """Set up the basic information needed by the parser

        Args:
            file_path:    Path to file that will be read.
            num_degrees:  Number of degrees of coefficients to read.
            encoding:     Encoding of file that will be read.
        """
        super().__init__(file_path, encoding=encoding)
        self.num_degrees = num_degrees
        self.raw = None

    def read_data(self):
        with open(self.file_path, mode="rt", encoding=self.file_encoding) as fid:
            self._read_header(fid)
            self._read_coeffs(fid)

        self._organize_data()

    def _read_header(self, fid):
        header_fields = {h.name: h for h in HEADER_FIELDS}
        for line in fid:
            if line.startswith("end_of_head"):
                break
            fields = line.strip().split()
            if fields and fields[0] in header_fields:
                header_def = header_fields[fields[0]]
                self.meta[header_def.name] = header_def.converter(fields[1])
        else:
            raise  # File ended without end_of_head

    def _read_coeffs(self, fid):
        if self.num_degrees is None:
            self.num_degrees = self.meta["max_degree"]

        # Skip orders below 2
        fid = itertools.dropwhile(lambda line: int(line.split()[1]) < 2, fid)

        # Read coefficient lines
        num_coeffs = int((self.num_degrees + 1) * (self.num_degrees + 2) / 2) - 3  # Skipping orders 00, 10, 11
        coeff_lines = itertools.islice(fid, num_coeffs)
        self.raw = np.genfromtxt(
            coeff_lines,
            names=("key", "degree", "order", "C", "S", "sigma_C", "sigma_S"),
            dtype=("U4", "i8", "i8", "f8", "f8", "f8", "f8"),
        )

    def _organize_data(self):
        self.data["C"] = np.zeros((self.num_degrees + 1, self.num_degrees + 1))
        self.data["S"] = np.zeros((self.num_degrees + 1, self.num_degrees + 1))
        self.data["C"][0, 0] = 1  # C[0, 0] is always 1

        # Add data to matrices
        self.data["C"][self.raw["degree"], self.raw["order"]] = self.raw["C"]
        self.data["S"][self.raw["degree"], self.raw["order"]] = self.raw["S"]
