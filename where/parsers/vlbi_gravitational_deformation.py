"""A parser for reading the gravitational deformation file
"""

# Standard library imports
from datetime import datetime

# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where.parsers._parser import Parser


@plugins.register
class GravitationalDeformationParser(Parser):
    """A parser for reading the gravitational deformation file
    """

    def read_data(self):
        with open(self.file_path, mode="rt") as fid:
            fid.readline()  # skip first line which only contains file version information
            current_station = ""

            for line in fid:
                if line.startswith("#") or not line:
                    continue

                line = line.split()
                if len(line) == 3:
                    # block header line
                    if line[0] == "EPOCH":
                        self.data[current_station]["start"] = datetime.strptime(line[1], "%Y%M%d")
                        self.data[current_station]["end"] = datetime.strptime(line[1], "%Y%M%d")
                    else:
                        current_station = line[0].strip().replace(" ", "_")
                        sta_dict = self.data.setdefault(current_station, dict())
                        sta_dict["elevation"] = list()
                        sta_dict["delay"] = list()
                        sta_dict["scale"] = float(line[2])
                        sta_dict["start"] = datetime.min
                        sta_dict["end"] = datetime.max
                else:
                    # data line
                    self.data[current_station]["elevation"].append(float(line[0]))  # Degrees
                    self.data[current_station]["delay"].append(float(line[1]))

    def calculate_data(self):
        for sta_dict in self.data.values():
            sta_dict["delay"] = np.array(sta_dict["delay"]) * sta_dict["scale"] * Unit.mm2m
            del sta_dict["scale"]
            sta_dict["elevation"] = np.deg2rad(sta_dict["elevation"])
