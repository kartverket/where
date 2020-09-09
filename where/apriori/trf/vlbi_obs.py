"""A TRF Factory for creating sites with positions from VLBI observation files

Description:
------------

Reads station positions and velocities from the VLBI observation files in either NGS or vgosdb format. No velocities
are available in the observation files, so the same position is used for all time epochs.
"""

# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import ellipsoid

# Where imports
from where.data.position import Position
from where import apriori
from where.apriori.trf import TrfFactory
from where.lib import config
from where.lib import log
from where import parsers


@plugins.register
class VlbiObsTrf(TrfFactory):
    """A factory for using positions from VLBI observation files
    """

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.version> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        fmt = config.tech.get("obs_format", self.version).str
        try:
            return getattr(self, "_read_data_{}".format(fmt))()
        except AttributeError:
            log.fatal(f"Format {fmt!r} is unknown for reference frame {self}")

    def _read_data_vgosdb(self):
        """Read data from vgosdb observation files

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        station_codes = apriori.get("vlbi_station_codes")
        data = parsers.parse_key("vlbi_obs_stations_vgosdb").as_dict()
        try:
            stations = [s.strip().replace(" ", "_") for s in data["AprioriStationList"]]
            xyz = data["AprioriStationXYZ"]
        except KeyError:
            return {}

        return {
            (station_codes[n]["cdp"] if n in station_codes else "key{}".format(i)): (
                dict(name=n, pos=p, **station_codes[n]) if n in station_codes else dict(name=n, pos=p, real=False)
            )
            for i, (n, p) in enumerate(zip(stations, xyz))
        }

    def _read_data_ngs(self):
        """Read data from NGS observation files

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        station_codes = apriori.get("vlbi_station_codes")
        pos_ngs = parsers.parse_key("vlbi_obs_ngs", parser_name="vlbi_ngs_sites").as_dict()

        # Match NGS names with station codes to get CDP codes as keys
        keys = {n: station_codes[n]["cdp"] for n in pos_ngs.keys() if n in station_codes}

        # Add info from station codes file to each site
        return {c: dict(name=n, pos=pos_ngs[n], **station_codes[n]) for n, c in keys.items()}

    def _calculate_pos_trs(self, site):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            site (String):  Key specifying which site to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """

        pos = self.data[site].pop("pos")[None, :].repeat(self.time.size, axis=0)
        ell = ellipsoid.get(config.tech.reference_ellipsoid.str.upper())
        pos_trs = Position(pos, system="trs", ellipsoid=ell, time=self.time)

        return np.squeeze(pos_trs)
