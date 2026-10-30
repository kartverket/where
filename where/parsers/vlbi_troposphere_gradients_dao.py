"""A parser for reading vlbi dao troposphere gradients

Description:
------------

Reads troposphere gradients precomputed by Nasa data assimilation office for vlbi stations

"""
# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class VlbiDaoTroposphereGradientsParser(LineParser):
    """A parser for reading vlbi dao troposphere gradients
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        # No header or comments in the file
        return dict(
            names=[
                "station",
                "gn",
                "gn_rms",
                "ge",
                "ge_rms"
            ],
            dtype=[
                "U8",
                "f8",
                "f8",
                "f8",
                "f8",
                ],
            comments="*",
            skip_header=1,
        )

    def structure_data(self):
        self.data = {
            item["station"]: dict(
                gn=item["gn"],
                ge=item["ge"],
                gn_rms=item["gn_rms"],
                ge_rms=item["ge_rms"],
            )
            for item in self._array
        }
