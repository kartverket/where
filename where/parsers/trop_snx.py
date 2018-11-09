"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads epoch of discontinuity for the position and velocity model in ITRF.

The file is using the SINEX format, but the SOLUTION/DISCONTINUITY block is not defined in the official format
description.

"""

# Standard library imports
from datetime import datetime

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_sinex import SinexParser, SinexBlock, SinexField, parsing_factory


@plugins.register
class TropSnxParser(SinexParser):
    """A parser for reading data from TRO files in SNX format
    """

    def setup_parser(self):
        #return {self.file_reference,self.trop_description,self.trop_sta_coordinates,self.trop_solution}
        return {self.trop_sta_coordinates,self.trop_solution}

    @property
    def trop_description(self):
        """Custom made block for ITRF to mark the epoch of discontinuities in the position and velocity of stations

        Example:
             1515  A    1 R 00:000:00000 92:180:43054 P - EQ M7.3 - Southern California
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
            
             SOLUTION_FIELDS_1             TROTOT STDDEV TGNTOT STDDEV TGETOT STDDEV
        """
        return SinexBlock(
            marker="TROP/DESCRIPTION",
            fields=(
                SinexField("keyword", 1, "U30"),
                SinexField("values", 31, "U49"),
            ),
            parser=self.parse_trop_description,
        )
        
    parse_trop_description=parsing_factory() 
    
    @property
    def trop_sta_coordinates(self):
        """Custom made block for ITRF to mark the epoch of discontinuities in the position and velocity of stations

        Example:
             1515  A    1 R 00:000:00000 92:180:43054 P - EQ M7.3 - Southern California
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
            
             ABMF  A    1 P  2919785.771 -5383744.975  1774604.831 IGS14  COD
        """
   
        return SinexBlock(
            marker="TROP/STA_COORDINATES",
            fields=(
                SinexField("site_code", 1, "U4"),
                SinexField("point_code", 6, "U2"),
                SinexField("soln", 9, "U4"),
                SinexField("obs_code", 14, "U1"),
                SinexField("sta_x", 16, "f8"),
                SinexField("sta_y", 29, "f8"),
                SinexField("sta_z", 42, "f8"),
                SinexField("trf_system", 55, "U6"),
                SinexField("remark", 62, "U5"),
            ),
            parser=self.parse_trop_sta_coordinates,
        )

    def parse_trop_sta_coordinates(self, data):
   
        for d in data:
            site_key = d["site_code"]
            self.data.setdefault(site_key, dict())
            self.data[site_key].update(dict(
                pos=np.array((d["sta_x"],d["sta_y"],d["sta_z"])),trf_system=d["trf_system"]
                ))
    
    @property
    def trop_solution(self):
        """Custom made block for ITRF to mark the epoch of discontinuities in the position and velocity of stations

        Example:
             1515  A    1 R 00:000:00000 92:180:43054 P - EQ M7.3 - Southern California
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
             ABMF 18:002:03600 2540.5    0.4  -0.437  0.039  -1.383  0.043
             
        """
        return SinexBlock(
            marker="TROP/SOLUTION",
            fields=(
                SinexField("site_code", 1, "U4"),
                SinexField("epoch", 6, "O","epoch"),
                SinexField("trotot", 19, "f8"),
                SinexField("stddev", 26, "f8"),
                SinexField("tgntot", 33, "f8"),
                SinexField("stddevn", 40, "f8"),
                SinexField("tgetot", 47, "f8"),
                SinexField("stddeve", 55, "f8"),
            ),
            parser=self.parse_trop_solution,
        )
        
    def parse_trop_solution(self, data):
        fields = ("epoch", "trotot", "stddev","tgntot","stddevn","tgetot","stddeve")
        sites = set()
   
        for d in data:
            site_key = d["site_code"]
            sites.add(site_key)
            self.data.setdefault(site_key, dict())
            for field in fields:
                self.data[site_key].setdefault(field, list()).append(d[field])

        for site in sites:
            for field in fields:
                self.data[site][field] = np.array(self.data[site][field])
              
