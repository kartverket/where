"""Get apriori data for eccentricities

Description:

Reads data for eccentricity files

"""
# Standard library imports
from collections import UserDict
from datetime import datetime, time

# Third party imports
import numpy as np

# Where imports
from where import parsers
from where.lib import plugins
from where.lib import log


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
        self._warned_missing = set()

    def _pick_data(self, data, rundate):
        for site_id, site_data in data.items():
            for key, value in site_data.items():
                if isinstance(key, tuple):
                    start = key[0]
                    end = key[1]
                    if datetime.combine(rundate, time.max) > start and datetime.combine(rundate, time.min) < end:
                        self.data.setdefault(site_id, {}).update(value)
                else:
                    self.data.setdefault(site_id, {}).update({key: value})

    def __missing__(self, site_id):
        """Handle missing keys

        Give a warning and return a consistent value when eccentricity vectors are missing from the file.

        Args:
            site_id (String):  site_id for the missing eccentricity vector

        Returns:
            Dict:  Dummy information including a zero eccentricity vector
        """
        if site_id not in self._warned_missing:
            log.warn(f"Missing eccentricity data for site id {site_id}. Vector set to zero.")
            self._warned_missing.add(site_id)
        return dict(vector=np.zeros(3), name="No name", coord_type="XYZ")
