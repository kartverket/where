"""A parser for reading data from ITRF files in SNX format

Description:
------------

Reads station positions and velocities from ITRF files in SNX format.

"""
# Standard library imports
from datetime import datetime

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_sinex import SinexParser


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
            self.data[site_key]["pos_vel"][soln] = dict(start=d["start_epoch"], end=d["end_epoch"], mean=d["mean_epoch"])

    def parse_solution_estimate(self, data):
        for d, soln in zip(data, data["soln"].astype("i8")):
            site_key = d["site_code"]
            self.data[site_key]["ref_epoch"] = d["ref_epoch"]
            self.data[site_key]["pos_vel"][soln][d["param_name"]] = d["estimate"]
