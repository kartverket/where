"""A parser for reading data from TRF files in SSC format

Description:
------------

Reads station positions and velocities from TRF files in SSC format. The velocity model is a simple linear offset
based on the reference epoch.

"""

# Standard library imports
from datetime import datetime, timedelta
import itertools
import re

# Midgard imports
from midgard.dev import plugins
from midgard.parsers._parser_chain import ParserDef, ChainParser


@plugins.register
class TrfSscParser(ChainParser):
    """A parser for reading data from ITRF files in SSC format
    """

    def setup_parser(self):
        # Read reference epoch from first line
        epoch_parser = ParserDef(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda _l, _ln: True,
            parser_def={True: {"parser": self.parse_epoch, "fields": {"ref_epoch": (0, 100)}}},
        )

        # Ignore header
        header_parser = ParserDef(
            end_marker=lambda _l, _ln, nextline: nextline[0:5].isnumeric(), label=None, parser_def=None
        )

        # Every pair of lines contains information about one station
        obs_parser = ParserDef(
            end_marker=lambda line, _ln, _n: line[30:31] == " ",
            label=lambda line, _ln: not line[30:31] == " ",
            parser_def={
                True: {
                    "parser": self.parse_position,
                    "fields": {
                        "site_num": (0, 5),
                        "antenna_num": (5, 9),
                        "name": (10, 26),
                        "tech": (26, 32),
                        "antenna_id": (32, 37),
                        "data": (37, 200),  # Column numbers for data are inconsistent
                    },
                },
                False: {"parser": self.parse_velocity, "fields": {"data": (37, 200)}},
            },
        )

        return itertools.chain([epoch_parser, header_parser], itertools.repeat(obs_parser))

    def parse_epoch(self, line, _):
        """Parse the epoch from the header line

        Looks for a number of the form yyyy.0, and stores it as the reference epoch.
        """
        epoch = re.search(r" (\d{4})\.0", line["ref_epoch"])
        if epoch:
            self.meta["ref_epoch"] = datetime.strptime(epoch.groups()[0], "%Y")

    def parse_position(self, line, cache):
        """Parse the position line of ITRF data

        This gives the position (x,y,z) of the station. Converting position float.

        Args:
            line (Dict):  The fields of a line.
            cache (Dict): Dict that persists information.
        """
        data_fields = ("STAX", "STAY", "STAZ", "sigma_x", "sigma_y", "sigma_z", "soln", "start", "end")
        data_values = line.pop("data")
        line.update({k: v for k, v in itertools.zip_longest(data_fields, data_values.split())})
        line["ref_epoch"] = self.meta["ref_epoch"]
        cache["antenna_id"] = line["antenna_id"]
        line["soln"] = int(line["soln"]) if line["soln"] else 1
        cache["soln"] = line["soln"]
        pos_vel = dict()
        pos_vel.update({k: float(line.pop(k)) for k in list(line.keys()) if k.startswith("STA")})

        start = line.pop("start")
        if start and start[3:6] != "000":
            pos_vel["start"] = datetime.strptime(start[0:6], "%y:%j") + timedelta(seconds=int(start[7:]))
        else:
            pos_vel["start"] = datetime.min

        end = line.pop("end")
        if end and end[3:6] != "000":
            pos_vel["end"] = datetime.strptime(end[0:6], "%y:%j") + timedelta(seconds=int(end[7:]))
        else:
            pos_vel["end"] = datetime.max

        self.data.setdefault(cache["antenna_id"], dict())
        self.data[cache["antenna_id"]].update(line)
        self.data[cache["antenna_id"]].setdefault("pos_vel", dict())
        self.data[cache["antenna_id"]]["pos_vel"][line["soln"]] = pos_vel

    def parse_velocity(self, line, cache):
        """Parsing the velocity line of ITRF data

        This is given on the line below the line with the position.  Assume that the tech and antenna_id are the
        same as on the line above, so we don't parse this.

        Args:
            line (Dict):  The fields of a line.
            cache (Dict): Dict that persists information.
        """
        data_fields = ("VELX", "VELY", "VELZ")
        data = {k: float(v) for k, v in zip(data_fields, line["data"].split())}
        self.data[cache["antenna_id"]]["pos_vel"][cache["soln"]].update(data)
