"""A parser for reading post-seismic deformation data from ITRF files in SNX format

Description:
------------

Reads post-seismic deformation model paramters from ITRF files in SNX format.

"""

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_sinex import SinexParser


@plugins.register
class TrfSnxPsdParser(SinexParser):
    """A parser for reading data from ITRF files in SNX format
    """

    def setup_parser(self):
        return (self.site_id, self.solution_estimate)

    def parse_site_id(self, data):
        for d in data:
            site_key = d["site_code"]
            self.data.setdefault(site_key, dict())
            self.data[site_key] = dict(
                antenna_id=d["site_code"], marker=d["marker"], domes=d["domes"], name=d["description"]
            )

    def parse_solution_estimate(self, data):
        for d in data:
            site_key = d["site_code"]
            self.data[site_key].setdefault("psd", dict()).setdefault(d["ref_epoch"], dict(epoch=d["ref_epoch"]))
            self.data[site_key]["psd"][d["ref_epoch"]].setdefault(d["param_name"], list()).append(d["estimate"])
