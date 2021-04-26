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

# Where imports
from where.lib import log


class HeaderField(NamedTuple):
    name: str
    converter: Callable


HEADER_FIELDS = (
    HeaderField("product_type", str),
    HeaderField("modelname", str),
    HeaderField("earth_gravity_constant", float),
    HeaderField("radius", float),
    HeaderField("max_degree", int),
    HeaderField("errors", str),
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

    def read_data(self):
        with open(self.file_path, mode="rt", encoding=self.file_encoding) as fid:
            self._read_header(fid)
            self._read_coeffs(fid)

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

        # Skip high orders
        fid = itertools.filterfalse(
            lambda line: (int(line.split()[1]) > self.num_degrees) or (int(line.split()[2]) > self.num_degrees), fid
        )

        # The EGM2008 file has some numbers in "fortran" format. Converter to correct for this:
        converter = lambda s: 1.0 if s == "1.0d0" else (0.0 if s == "0.0d0" else float(s))

        # Read coefficient lines
        for k, coeff_lines in itertools.groupby(sorted(fid), lambda line: line.split()[0]):

            if k == "gfc":
                self.data[k] = np.genfromtxt(
                    coeff_lines,
                    names=("key", "degree", "order", "C", "S", "sigma_C", "sigma_S"),
                    dtype=("U4", "i8", "i8", "U30", "U30", "U30", "U30"),
                    encoding="UTF-8",
                    converters={3: converter, 4: converter, 5: converter, 6: converter},
                )
            elif k == "trnd":
                self.data[k] = np.genfromtxt(
                    coeff_lines,
                    names=("key", "degree", "order", "C", "S", "sigma_C", "sigma_S", "t0", "t1"),
                    dtype=("U4", "i8", "i8", "f8", "f8", "f8", "f8", "U13", "U13"),
                )
            elif k == "gfct" and self.meta["errors"] in ("formal", "calibrated"):
                self.data[k] = np.genfromtxt(
                    coeff_lines,
                    names=("key", "degree", "order", "C", "S", "sigma_C", "sigma_S", "t0", "t1"),
                    dtype=("U4", "i8", "i8", "f8", "f8", "f8", "f8", "U13", "U13"),
                )
            elif k in ("asin", "acos") and self.meta["errors"] in ("formal", "calibrated"):
                self.data[k] = np.genfromtxt(
                    coeff_lines,
                    names=("key", "degree", "order", "C", "S", "sigma_C", "sigma_S", "t0", "t1", "period"),
                    dtype=("U4", "i8", "i8", "f8", "f8", "f8", "f8", "U13", "U13", "f8"),
                )
            else:
                log.warn("Unknown gravity coefficient type")
