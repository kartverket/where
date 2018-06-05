"""A parser for reading radio source coordinates from VASCC apriori crf

Description:
------------

Reads radio source coordinates from VASCC (VLBI Software Analysis Comparison Campaign) apriori file.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# External library imports
import numpy as np

# Where imports
from where.parsers._parser_line import LineParser
from where.lib import plugins
from where.lib import files


@plugins.register
class VasccCrfParser(LineParser):
    """A parser for reading source coordinates from ICRF files
    """

    def setup_parser(self):
        return dict(dtype=("U8", "U4", "f8", "f8", "f8", "f8", "f8", "f8"), skip_header=2)

    def structure_data(self):
        ref_epoch = "2005.0"
        self.data = {
            cdp: {"antenna_id": cdp, "ref_epoch": ref_epoch, "site_name": name, "pos": [x, y, z], "vel": [vx, vy, vz]}
            for name, cdp, x, y, z, vx, vy, vz in self._array
        }
