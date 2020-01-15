"""A parser for reading gLAB output files

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='glab_output', file_path='glab_output.txt')
    data = p.as_dict()

Description:
------------


"""
# Standard library imports
from datetime import datetime, timedelta
from typing import Any, Dict, List, Union

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins
from midgard.parsers.glab_output import GlabOutputParser


@plugins.register
class GlabOutputParser(GlabOutputParser):
    """A parser for reading gLAB output files

    The keys of the **data** dictionary are defined depending, which kind of gLAB output file is read. The values of 
    the **data** dictionary are represented by the gLAB colum values.

    Following **meta**-data are available after reading of gLAB files:

    | Key                  | Description                                                                          |
    |----------------------|--------------------------------------------------------------------------------------|
    | \__data_path__       | File path                                                                            |
    | \__parser_name__     | Parser name                                                                          |
    """

    def write_to_dataset(self, dset: "Dataset") -> "Dataset":
        """Return the parsed data as a Dataset

        Args:
            dset: A dataset containing the data.

        Returns:
            A dataset containing the data.
        """

        # Initialize dataset
        if not self.data:
            log.warn("No data in {self.file_path}.")
            return dset
        dset.num_obs = len(self.data["year"])

        # Add time
        epochs = list()
        for year, doy, seconds in zip(self.data["year"], self.data["doy"], self.data["seconds"]):
            epochs.append(datetime.strptime("{:.0f} {:.0f}".format(year, doy), "%Y %j") + timedelta(seconds=seconds))

        dset.add_time("time", val=epochs, scale="gps", format="datetime", write_level="operational")

        # Add system field
        if "system" in self.data.keys():
            systems = []
            for system in self.data["system"]:
                systems.append(enums.gnss_name_to_id[system.lower()].value)

            dset.add_text("system", val=systems)

        # Add system field
        if "satellite" in self.data.keys():
            satellites = []
            for system, satellite in zip(dset.system, self.data["satellite"]):
                satellites.append(system + str(satellite).zfill(2))

            dset.add_text("satellite", val=satellites)

        # Add text and float fields
        fields = set(self.data.keys()) - {"year", "doy", "seconds", "system", "satellite"}
        for field in fields:
            if self.data[field].dtype.kind in {"U", "S"}:  # Check if numpy type is string
                dset.add_text(field, val=self.data[field])
                continue

            dset.add_float(field, val=self.data[field])

        return dset
