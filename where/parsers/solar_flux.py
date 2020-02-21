"""A parser for solar flux from file

Example:
    from where.parsers import SolarFluxParser
    parser = SolarFluxParser(file_key)
    parser.process_data()

Description:

    Reads data from reference below, and stores in table.

References:
     http://www.ngdc.noaa.gov/stp/space-weather/solar-data/solar-features/solar-radio/noontime-flux/penticton/penticton_observed/tables/

@todo All these files from the link above are merged to F10.7CM
           (leftovers from old GEOSAT). Do something about this?

"""

# Standard library imports
from datetime import datetime
import itertools

# Midgard imports
from midgard.dev import plugins
from midgard.math.constant import constant
from midgard.parsers._parser_chain import ChainParser, ParserDef


@plugins.register
class SolarFluxParser(ChainParser):
    """A parser for reading Solar Flux data
    """

    def setup_parser(self):
        flux_parser = ParserDef(
            end_marker=lambda _l, _ln, nextline: nextline[0:1].isdigit(),
            label=lambda line, _ln: line[0:1].isdigit(),
            parser_def={
                True: {"parser": self.parse_header, "fields": {"year": (0, 4)}},
                False: {
                    "parser": self.parse_flux,
                    "fields": {
                        "day": (3, 5),
                        "01": (6, 11),
                        "02": (12, 17),
                        "03": (18, 23),
                        "04": (24, 29),
                        "05": (30, 35),
                        "06": (36, 41),
                        "07": (42, 47),
                        "08": (48, 53),
                        "09": (54, 59),
                        "10": (60, 65),
                        "11": (66, 71),
                        "12": (72, 77),
                    },
                },
            },
        )

        return itertools.repeat(flux_parser)

    def parse_header(self, line, cache):
        cache["year"] = int(line.pop("year"))

    def parse_flux(self, line, cache):
        if not line["day"]:
            return
        day = int(line.pop("day"))
        for month in line:
            if not line[month]:
                continue
            time = None
            try:
                time = datetime(cache["year"], int(month), day)
                self.data[time] = int(line[month])
            except ValueError:  # Todo: Can we use __missing__ instead? Makes it possible to add warning
                if time is not None:
                    self.data[time] = constant.get("S", source="book")
