"""A parser for reading Terrapos residual file

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='terrapos_residual', file_path='PPP-residuals.txt')
    data = p.as_dict()

Description:
------------

Reads data from files in Terrapos residual format.

"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_line import LineParser


@plugins.register
class TerraposResidualParser(LineParser):
    """A parser for reading Terrapos residual file

    Following parameters are available after reading Terrapos residual file:

        ====================  =================================================================
         Parameter             Description
        ====================  =================================================================
         azimuth               Azimuth of satellites in [deg]
         elevation             Elevation of satellites in [deg]
         gpsweek               GPS week
         gpssec                Seconds of GPS week
         residual_code         Code (pseudorange) residuals in [m]
         residual_doppler      Doppler residuals in [m]
         residual_phase        Carrier-phase residuals in [m]
         satellite             Satellite PRN number together with GNSS identifier (e.g. G07)
         system                GNSS identifier
        ====================  =================================================================
    """

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """

        # Sat Week    ToW (s)   Azi  Elev       code_res      phase_res    doppler_res
        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+-
        # G01 1972 518400.000  32.0   8.0         -0.560          0.038            NaN
        # G02 1972 518400.000  25.0   6.2          2.430          0.078            NaN
        return dict(
            names=(
                "satellite",
                "gpsweek",
                "gpssec",
                "azimuth",
                "elevation",
                "residual_code",
                "residual_phase",
                "residual_doppler",
            ),
            delimiter=(3, 5, 11, 6, 6, 15, 15, 15),
            dtype=("U3", "f8", "f8", "f8", "f8", "f8", "f8", "f8"),
        )

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return [self._add_system_field]

    def _add_system_field(self):
        """Add system parameter to data
        """
        self.data["system"] = self.data["satellite"].astype("U1")
