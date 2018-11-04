"""A parser for reading radio source coordinates from VASCC apriori crf

Description:
------------

Reads radio source coordinates from VASCC (VLBI Software Analysis Comparison Campaign) apriori file.

"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_line import LineParser


@plugins.register
class VasccCrfParser(LineParser):
    """A parser for reading source coordinates from ICRF files
    """

    def setup_parser(self):
        return dict(usecols=(0, 3, 4), dtype="U8, f8, f8", skip_header=1)

    def structure_data(self):
        self.data = {
            name: {
                "ra": ra,
                "dec": dec,
                "special": False,
                "undefined": True,
                "non_vcs": False,
                "vcs": False,
                "defining": False,
            }
            for name, ra, dec in self._array
        }
