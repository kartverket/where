"""Get apriori data for eccentricities

Description:

Reads data for eccentricity files

"""
# Standard library imports
from collections import UserDict
from datetime import datetime, time

# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import parsers


@plugins.register
def get_eccentricity(rundate):
    """Get Eccentricities for a given date

    Args:
        rundate: The run date of the data analysis.

    Returns:
        A dictionary of eccentricities.
    """

    data = parsers.parse_key(file_key="eccentricity").as_dict()
    return Eccentricity(data, rundate)


class Eccentricity(UserDict):
    def __init__(self, data, rundate):
        super().__init__()
        self._pick_data(data, rundate)

    def _pick_data(self, data, rundate):
        for identifier, site_data in data.items():
            for key, value in site_data.items():
                if isinstance(key, tuple):
                    start = key[0]
                    end = key[1]
                    if datetime.combine(rundate, time.max) > start and datetime.combine(rundate, time.min) < end:
                        self.data.setdefault(identifier, {}).update(value)
                else:
                    self.data.setdefault(identifier, {}).update({key: value})

    def __missing__(self, identifier):
        """Handle missing keys

        Args:
            identifier (String):  identifier for the site

        Returns:
            Dict:  Dummy information including a zero eccentricity vector
        """
        return dict(vector=np.zeros(3), name="No name", coord_type="XYZ")
