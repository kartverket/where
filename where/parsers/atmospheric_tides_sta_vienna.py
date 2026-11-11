"""A parser for reading atmospheric tide coefficients

Description:
------------

Reads atmospheric tide coefficients from Vienna. This model is expressed in a center of mass (cm) frame.

Station specific coefficients for the S1, S2 and S3 atmospheric tides.

Column description tidal APL data: 
(1)      station name
(2)      S1 cosine amplitude for radial displacement [mm]
(3)      S1 sine amplitude for radial displacement [mm]
(4)      S2 cosine amplitude for radial displacement [mm]
(5)      S2 sine amplitude for radial displacement [mm]
(6)      S3 cosine amplitude for radial displacement [mm]
(7)      S3 sine amplitude for radial displacement [mm]
(8)      S1 cosine amplitude for eastward displacement [mm]
(9)      S1 sine amplitude for eastward displacement [mm]
(10)     S2 cosine amplitude for eastward displacement [mm]
(11)     S2 sine amplitude for eastward displacement [mm]
(12)     S3 cosine amplitude for eastward displacement [mm]
(13)     S3 sine amplitude for eastward displacement [mm]
(14)     S1 cosine amplitude for northward displacement [mm]
(15)     S1 sine amplitude for northward displacement [mm]
(16)     S2 cosine amplitude for northward displacement [mm]
(17)     S2 sine amplitude for northward displacement [mm]
(18)     S3 cosine amplitude for northward displacement [mm]
(19)     S3 sine amplitude for northward displacement [mm]


"""
# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_line import LineParser


@plugins.register
class ViennaStaAtmosphericTidesParser(LineParser):
    """A parser for reading stationwise coefficents for atmospheric tides
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
                "A_d1_u",
                "B_d1_u",
                "A_d2_u",
                "B_d2_u",
                "A_d3_u",
                "B_d3_u",
                "A_d1_e",
                "B_d1_e",
                "A_d2_e",
                "B_d2_e",
                "A_d3_e",
                "B_d3_e",
                "A_d1_n",
                "B_d1_n",
                "A_d2_n",
                "B_d2_n",
                "A_d3_n",
                "B_d3_n",
            ],
            dtype = ["U8"] + ["f8"]*18,
            comments="#",
        )

    def structure_data(self):
        self.data = {
            item["station"]: {k:item[k] for k in item.dtype.names[1:]}
            for item in self._array
        }
