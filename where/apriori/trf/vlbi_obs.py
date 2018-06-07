"""A TRF Factory for creating sites with positions from VLBI observation files

Description:
------------

Reads station positions and velocities from the VLBI observation files in either NGS or vgosdb format. No velocities
are available in the observation files, so the same position is used for all time epochs.



"""

import netCDF4

# Where imports
from where import apriori
from where.apriori import trf
from where.lib import config
from where.lib import log
from where.lib import plugins
from where import parsers
from where.lib import files
from where.lib import util


@plugins.register
class VlbiObsTrf(trf.TrfFactory):
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
            log.fatal("Format '{}' is unknown for reference frame {}", fmt, self)

    def _read_data_vgosdb(self):
        """Read data from vgosdb observation files

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        # TODO create parser
        stations_file = files.path(file_key="vlbi_obs_vgosdb") / "Apriori" / "Station.nc"
        d = netCDF4.Dataset(stations_file)
        stations = d["AprioriStationList"][:]
        xyz = d["AprioriStationXYZ"][:]
        station_codes = apriori.get("vlbi_station_codes")
        stations = [str(n, "utf-8").strip() for n in stations]
        return {
            station_codes[n]["cdp"]: dict(name=n, pos=p, **station_codes[n])
            for n, p in zip(stations, xyz)
            if n in station_codes
        }

    def _read_data_ngs(self):
        """Read data from NGS observation files

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        station_codes = apriori.get("vlbi_station_codes")
        pos_ngs = parsers.parse_key("vlbi_obs_ngs", parser="vlbi_ngs_sites").as_dict()

        # Match NGS names with station codes to get CDP codes as keys
        keys = {n: station_codes[n]["cdp"] for n in pos_ngs.keys() if n in station_codes}

        # Add info from station codes file to each site
        return {c: dict(name=n, pos=pos_ngs[n], **station_codes[n]) for n, c in keys.items()}

    def _calculate_pos_itrs(self, site):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            site (String):  Key specifying which site to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        return self.data[site]["pos"][None, :].repeat(len(self.time), axis=0)
