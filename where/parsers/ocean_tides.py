"""A parser for reading ocean tidal loading coefficients

Example:
    from where.parsers import OceanTidesParser
    parser = OceanTidesParser()
    parser.process_data()

Description:

Reads data from files in the BLQ file format provided by the ocean tide
loading service at http://holt.oso.chalmers.se/loading.

Note: The format is not defined explicitly, but is illustrated by the
example http://holt.oso.chalmers.se/loading/example_blq.html.

The following assumption is made by this parser:
*) A data block is identified with a 5 digit domes number

"""

# Standard library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser import Parser


@plugins.register
class OceanTidesParser(Parser):
    """A parser for reading ocean tidal loading coefficients from BLQ-files
    """

    def read_data(self):
        with open(self.file_path, mode="rt") as fid:
            self._parse_file(fid)

    def _parse_file(self, fid):
        for line in fid:
            if line.startswith("$$") or not line:
                continue

            site_id = line.split()[0]
            self._parse_block(fid, site_id)

    def _parse_block(self, fid, site_id):
        self.data[site_id] = dict()
        self._parse_data(fid, site_id, "amplitudes")
        self._parse_data(fid, site_id, "phases")

    def _parse_data(self, fid, site_id, label):
        up = self._read_line(fid)
        west = self._read_line(fid)
        south = self._read_line(fid)
        self.data[site_id][label] = np.hstack((up, west, south)).reshape((3, -1))

    def _read_line(self, fid):
        for line in fid:
            if not line.startswith("$$"):
                break
        line = line.strip().split()
        return [float(f) for f in line]
