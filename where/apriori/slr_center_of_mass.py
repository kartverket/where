"""Description:

    Reads station-dependent information from file.  This is needed in order to compute center of mass corrections for
    the satellites lageos-1, lageos-2, etalon-1, etalon-2 and ajisai.

References:

    http://ilrs.dgfi.tum.de

"""

# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log
from where import parsers
from where.lib.unit import Unit

DEFAULT_COM = {"lageos": 0.251, "etalon": 0.65, "lares": 0.13}


@plugins.register
def get_center_of_mass(sat_name):
    """Read station-dependent center of mass corrections from file

    Args:
        sat_name (String):  Name of satellite.

    Returns:
        Dict:  Center of mass corrections per station for the given satellite.
    """
    satellite = sat_name.rstrip("12")
    parser = parsers.parse_key("slr_center_of_mass", file_vars=dict(satellite=satellite))
    return SlrCenterOfMass(parser.as_dict(), satellite)


class SlrCenterOfMass:
    def __init__(self, data, satellite):
        self.data = data
        self.satellite = satellite

    def get(self, stations, time):
        """Get center of mass corrections for the given stations and epochs

        Args:
            station:  Array of stations.
            time:     Array of time epochs.

        Returns:
            Array of center of mass corrections in meters.
        """
        com = np.ones(stations.shape) * DEFAULT_COM[self.satellite]
        for station in np.unique(stations):
            if station not in self.data:
                log.warn(
                    f"Missing center of mass data for CDP '{station}'. "
                    f"Using default for {self.satellite} ({DEFAULT_COM[self.satellite]} meters)."
                )
                continue
            idx = stations == station
            epochs = time[idx].utc.datetime
            station_com = np.full(sum(idx), np.nan, dtype=float)
            for interval in self.data[station]:
                epoch_idx = (interval["start"] <= epochs) & (epochs <= interval["end"])
                station_com[epoch_idx] = interval["lcm"] * Unit.millimeter2meter
            if np.any(np.isnan(station_com)):
                first_epoch = epochs[np.where(np.isnan(station_com))[0][0]]
                log.warn(
                    f"Missing center of mass data for CDP '{station}' for {first_epoch:%Y-%m-%d}. "
                    f"Using default for {self.satellite} ({DEFAULT_COM[self.satellite]} meters)."
                )
                station_com[np.isnan(station_com)] = DEFAULT_COM[self.satellite]
            com[idx] = station_com
        return com
