"""A parser for reading VLBI data from NGS files

Description:
------------

Reads data from files in the NGS file format as defined in http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt
(revision date June 11, 2007) [1].

References:
-----------

[1] NGS file format.
    http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt

[2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010).
    http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html



"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import constant
from where.lib import log
from where.parsers._parser_chain import ParserDef, ChainParser
from where.lib.unit import unit


@plugins.register
class VlbiNgs2Parser(ChainParser):
    """A parser for reading VLBI data from NGS files
    """

    def setup_parser(self):
        # Ignore the header, first two lines
        header_parser = ParserDef(end_marker=lambda _l, line_num, _n: line_num == 2, label=None, parser_def=None)

        station_parser = ParserDef(end_marker=lambda line, _ln, _n: line == "$END", label=None, parser_def=None)
        source_parser = ParserDef(end_marker=lambda line, _ln, _n: line == "$END", label=None, parser_def=None)
        param_parser = ParserDef(end_marker=lambda line, _ln, _n: line == "$END", label=None, parser_def=None)

        # Observations are listed on 9 lines
        obs_parser = ParserDef(
            end_marker=lambda _l, _ln, next_line: next_line[78:80] == "01",
            # end_callback=self.copy_cache_to_obs,
            label=lambda line, _ln: line[78:80],
            parser_def={
                "01": {
                    "parser": self.parse_obs_meta,
                    "fields": {
                        "station_1": (0, 8),
                        "station_2": (10, 18),
                        "source": (20, 28),
                        "year": (29, 33),
                        "month": (34, 36),
                        "day": (37, 39),
                        "hour": (40, 42),
                        "minute": (43, 45),
                        "seconds": (46, 60),
                    },
                },
                "02": {
                    "parser": self.parse_obs(unit_in="nanoseconds", except_fields=("data_quality",)),
                    "fields": {
                        "observed_delay": (0, 20),
                        "observed_delay_ferr": (20, 30),
                        #                                  'observed_delay_rate':      (30, 50),
                        #                                  'observed_delay_rate_ferr': (50, 60),
                        "data_quality": (60, 62),
                        #                                  'flag_delay_type':          (63, 65),
                        #                                  'flag_delay_rate_type':     (66, 68),
                    },
                },
                #                 "03": {
                #                     "parser": self.parse_obs(),
                #                     "fields": {
                #                         "correlation_coeff": (0, 10),
                #                         "correlation_coeff_ferr": (10, 20),
                #                         "fringe_amplitude": (20, 30),
                #                         "fringe_amplitude_ferr": (30, 40),
                #                         "total_fringe_phase": (40, 60),
                #                         "total_fringe_phase_ferr": (60, 70),
                #                     },
                #                 },
                "05": {
                    "parser": self.parse_obs(unit_in="nanoseconds"),
                    "fields": {"cable_delay_1": (0, 10), "cable_delay_2": (10, 20)},
                },
                "06": {
                    "parser": self.parse_obs_missing,
                    "fields": {
                        "temperature_1": (0, 10),
                        "temperature_2": (10, 20),
                        "pressure_1": (20, 30),
                        "pressure_2": (30, 40),
                    },
                },
                "08": {
                    "parser": self.parse_obs(unit_in="nanoseconds", except_fields=("iono_quality",)),
                    "fields": {
                        "iono_delay": (0, 20),
                        "iono_delay_ferr": (20, 30),
                        #                                  'iono_delay_rate':           (30, 50),
                        #                                  'iono_delay_rate_ferr':      (50, 60),
                        "iono_quality": (61, 63),
                    },
                },
            },
        )

        return itertools.chain(
            [header_parser, station_parser, source_parser, param_parser], itertools.repeat(obs_parser)
        )

    def parse_obs_meta(self, line, cache):
        """Reads meta information like time stamp and station and source id

        Creates a new observation on the dataset with proper metainformation. Observation data are added to the dataset
        later by a calculator as it is more effective not to continuously resize numpy arrays.

        Args:
            line:  Input data from NGS file

        """
        line["seconds"] = "{:013.10f}".format(float(line["seconds"]))
        try:
            line["year"] = "{:4d}".format(int(line["year"]))
        except ValueError as e:
            # In some sessions the year field is '19 0' when it should be '2000'
            year = line["year"].split()
            if len(year) == 2:
                if int(year[1]) < 50:
                    line["year"] = "{:4d}".format(2000 + int(year[1]))
                else:
                    line["year"] = "{:4d}".format(1900 + int(year[1]))
            else:
                raise e

        obs = {
            "time": "{year:0>4}-{month:0>2}-{day:0>2}T{hour:0>2}:{minute:0>2}:{seconds}".format(**line),
            "station_1": line["station_1"],
            "station_2": line["station_2"],
            "source": line["source"],
        }

        for field, value in obs.items():
            self.data.setdefault(field, list()).append(value)

    def parse_obs(self, unit_in="meter", except_fields=()):
        """Read information about an observation

        Stores the information in the temporary cache-dict, which will be transfered to self.data when all information
        about this observation is parsed. If `unit_in` is given, all values will be converted to meter unless the field
        is listed in the `except_fields`-list.

        Args:
            unit_in (String):      Name of unit of values to be parsed.
            except_fields (Tuple): Names of fields where values should not be converted to meters.
        """
        # Find scale factor for converting to meter
        if unit_in == "meter":
            scale_factor = 1
        else:
            quantity = unit(unit_in)
            try:
                scale_factor = quantity.to("meter").magnitude
            except unit.DimensionalityError:
                # Try to convert between time and length by multiplying by the speed of light
                scale_factor = (quantity * constant.c * unit("meters per second")).to("meter").magnitude
        obs = {}

        # Define the function doing the actual parsing
        def parse_func(line, cache):
            for field in line:
                if field.startswith("flag_"):  # Flags are currently ignored
                    continue
                try:
                    obs[field] = float(line[field])
                except ValueError:
                    obs[field] = 0
                    log.debug("Could not convert {} to a number for {}. Value set to 0.0.", line[field], field)
                if field not in except_fields:
                    obs[field] *= scale_factor
            for field, value in obs.items():
                self.data.setdefault(field, list()).append(value)

        return parse_func

    def parse_obs_missing(self, line, cache):
        for key, value in line.items():
            if value.startswith("-999"):
                line[key] = "nan"
        self.parse_obs()(line, cache)
