"""A parser for reading SLR eccentricity vectors from file

Description:
------------

Reads the SLR eccentricity vector from file.

"""
# Standard library imports
from datetime import datetime

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_sinex import SinexParser, SinexBlock, SinexField


@plugins.register
class SlrEccentricityParser(SinexParser):
    """A parser for reading data from ITRF files in SNX format
    """

    max_line_width = 100

    def setup_parser(self):
        return (self.slr_site_id, self.slr_site_eccentricity)

    @property
    def slr_site_id(self):
        """General information for each site containing estimated parameters.

        Extra column added with 8 digit CDP SOD number

        Example:
            *CODE PT __DOMES__ T _STATION DESCRIPTION__ _LONGITUDE_ _LATITUDE__ HEIGHT_
             1515  A 40405S019 R DSS15    34-m HEF at G 243 06 46.0  35 25 18.6   973.2
                      1111111111222222222233333333334444444444555555555566666666667777777777
            01234567890123456789012345678901234567890123456789012345678901234567890123456789
        """
        return SinexBlock(
            marker="SITE/ID",
            fields=(
                SinexField("site_code", 1, "U4"),
                SinexField("point_code", 6, "U2"),
                SinexField("domes", 9, "U5"),
                SinexField("marker", 14, "U4"),
                SinexField("obs_code", 19, "U1"),
                SinexField("description", 21, "U22", "utf8"),
                SinexField("approx_lon", 44, "f8", "dms2deg"),
                SinexField("approx_lat", 56, "f8", "dms2deg"),
                SinexField("approx_height", 68, "f8"),
                SinexField("cdp_sod", 80, "U8"),
            ),
            parser=self.parse_site_id,
        )

    @property
    def slr_site_eccentricity(self):
        """List of antenna eccentricities

        Extra column added with 8 digit CDP SOD number

        Antenna eccentricities from the Marker to the Antenna Reference Point (ARP) or to the intersection of axis.
        """
        return SinexBlock(
            marker="SITE/ECCENTRICITY",
            fields=(
                SinexField("site_code", 1, "U4"),
                SinexField("point_code", 6, "U2"),
                SinexField("soln", 9, "U4"),
                SinexField("obs_code", 14, "U1"),
                SinexField("start_time", 16, "O", "epoch"),
                SinexField("end_time", 29, "O", "epoch"),
                SinexField("vector_type", 42, "U3"),
                SinexField("vector_1", 46, "f8"),
                SinexField("vector_2", 55, "f8"),
                SinexField("vector_3", 64, "f8"),
                SinexField("cdp_sod", 80, "U8"),
            ),
            parser=self.parse_site_eccentricity,
        )

    def parse_site_id(self, data):
        for d in data:
            site_key = d["site_code"]
            self.data.setdefault(site_key, dict())
            self.data[site_key] = dict(
                antenna_id=d["site_code"],
                marker=d["marker"],
                domes=d["domes"],
                name=d["description"],
                cdp_sod=d["cdp_sod"],
            )

    def parse_site_eccentricity(self, data):
        for d in data:
            start_time = datetime.min if d["start_time"] is None else d["start_time"]
            end_time = datetime.max if d["end_time"] is None else d["end_time"]
            key = (start_time, end_time)
            if d["vector_type"] == "UNE":
                # Convert UNE to ENU
                self.data[d["site_code"]].setdefault(key, {}).update(
                    dict(vector=(d["vector_3"], d["vector_2"], d["vector_1"]), coord_type="ENU")
                )
            else:
                self.data[d["site_code"]].setdefault(key, {}).update(
                    dict(vector=(d["vector_1"], d["vector_2"], d["vector_3"]), coord_type=d["vector_type"])
                )
