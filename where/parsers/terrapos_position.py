"""A parser for reading Terrapos position output file

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='terrapos_position', file_path='Gal_C1X_brdc_land_30sec_24hrs_FNAV-file.txt')
    data = p.as_dict()

Description:
------------

Reads data from files in Terrapos position output format.



"""
# Standard library imports
from typing import Dict

# Midgard imports
from midgard.parsers import LineParser
from midgard.dev import plugins


@plugins.register
class TerraposPositionParser(LineParser):
    """A parser for reading Terrapos position output file
    """

    def setup_parser(self) -> Dict=[str, str]:
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file. Following parameters are read:

        ====================  ======================================================================
         Parameter             Description
        ====================  ======================================================================
         gpsweek               GPS week
         gpssec                Seconds of GPS week
         head                  Head in [deg]
         height                Ellipsoidal height in [m]
         lat                   Latitude in [deg]
         lon                   Longitude in [deg]
         num_sat               Number of satellites
         pdop                  Position Dilution of Precision (PDOP)
         pitch                 Pitch in [deg]
         reliability_east      East position external reliability in [m] #TODO: Is that correct?
         reliability_height    Height position external reliability in [m] #TODO: Is that correct?
         reliability_north     North position external reliability in [m] #TODO: Is that correct?
         roll                  Roll in [deg]
         sigma_east            Standard deviation of East position in [m] #TODO: Is that correct?
         sigma_height          Standard deviation of Height position in [m] #TODO: Is that correct?
         sigma_north           Standard deviation of North position in [m] #TODO: Is that correct?
        ====================  ======================================================================
        """

        # Parser for Terrapos position output file

        # # Week       ToW (s)      Lat (deg)       Lon (deg)   Hght (m)    Roll(d)   Pitch(d)    Head(d)  sN (m)  sE (m)  sH (m)  sR (m)  sP (m) sHd (m)  #S    PDOP    rN (m)    rE (m)    rH (m)
        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+-
        #  1972  518430.000000   70.083249333    29.735357183    56.7994                                    5.300   5.895  10.821                           6     3.5   151.851   373.061   321.308
        #  1972  518460.000000   70.083251198    29.735329955    56.2342                                    5.303   5.925  10.778                           6     3.5   149.417   369.513   316.813
        return dict(
            comments="#",  # Remove comment lines starting with '#'
            names=(
                "gpsweek",
                "gpssec",
                "lat",
                "lon",
                "height",
                "roll",
                "pitch",
                "head",
                "sigma_north",
                "sigma_east",
                "sigma_height",
                "sigma_roll",
                "sigma_pitch",
                "sigma_head",
                "num_sat",
                "pdop",
                "reliability_north",
                "reliability_east",
                "reliability_height",
            ),
            delimiter=(5, 15, 15, 16, 11, 11, 11, 11, 8, 8, 8, 8, 8, 8, 4, 8, 10, 10, 10),
        )
