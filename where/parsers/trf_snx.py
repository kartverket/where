"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads station positions and velocities from ITRF files in SNX format.




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""
# Standard library imports
from datetime import datetime

# Where imports
from where.lib import plugins
from where.parsers._parser_sinex import SinexParser


@plugins.register
class TrfSnxParser(SinexParser):
    """A parser for reading data from ITRF files in SNX format
    """

    def setup_parser(self):
        return (self.site_id, self.solution_epochs, self.solution_estimate)

    def parse_site_id(self, data):
        for d in data:
            site_key = d["site_code"]
            self.data.setdefault(site_key, dict())
            self.data[site_key] = dict(
                antenna_id=d["site_code"], marker=d["marker"], domes=d["domes"], name=d["description"]
            )

    def parse_solution_epochs(self, data):
        for d, soln in zip(data, data["soln"].astype("i8")):
            site_key = d["site_code"]
            self.data[site_key].setdefault("pos_vel", dict())
            self.data[site_key]["pos_vel"][soln] = dict(start=datetime.min, end=datetime.max)

    def parse_solution_estimate(self, data):
        for d, soln in zip(data, data["soln"].astype("i8")):
            site_key = d["site_code"]
            self.data[site_key]["ref_epoch"] = d["ref_epoch"]
            self.data[site_key]["pos_vel"][soln][d["param_name"]] = d["estimate"]
