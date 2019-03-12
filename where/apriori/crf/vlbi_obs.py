"""A CRF Factory for creating sources with positions from VLBI observation files

Description:
------------

Reads source positions from the VLBI observation files in either NGS or vgosdb format.



"""

# External library imports
import numpy as np

# Where imports
from where import apriori
from where.apriori import crf
from where.lib import config
from where.lib import log
from where.lib import plugins
from where import parsers


@plugins.register
class VlbiObsCrf(crf.CrfFactory):
    """A factory for using source positions from VLBI observation files
    """

    def _read_data(self):
        """Read data needed by this Celestial Reference Frame for calculating positions of sources

        Delegates to _read_data_{<format>} to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each source defined in this celestial reference frame.
        """
        fmt = config.tech.obs_format.str
        try:
            return getattr(self, f"_read_data_{fmt}")()
        except AttributeError:
            log.fatal(f"Format {fmt!r} is unknown for reference frame {self}")

    def _read_data_vgosdb(self):
        """Read data from vgosdb observation files

        Returns:
            Dict:  Dictionary containing data about each source defined in this celestial reference frame.
        """
        source_names = apriori.get("vlbi_source_names")
        data = parsers.parse_key("vlbi_obs_sources_vgosdb").as_dict()
        try:
            sources = data["AprioriSourceList"]
            radec = data["AprioriSource2000RaDec"]
        except KeyError:
            return {}
        # Replace IVS name of source with official IERS name
        return {
            source_names[ivsname]["iers_name"] if ivsname in source_names else ivsname: dict(ra=coord[0], dec=coord[1])
            for ivsname, coord in zip(sources, radec)
        }

    def _read_data_ngs(self):
        """Read data from NGS observation files

        Returns:
            Dict:  Dictionary containing data about each source defined in this celestial reference frame.
        """
        data = parsers.parse_key("vlbi_obs_ngs", parser_name="vlbi_ngs_sources").as_dict()
        src_names = apriori.get("vlbi_source_names")

        # Replace IVS name of source with official IERS name
        return {
            src_names[ivsname]["iers_name"] if ivsname in src_names else ivsname: coords
            for ivsname, coords in data.items()
        }

    def _calculate_pos_crs(self, source):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            source (String):  Key specifying which source to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 2-vector.
        """
        return np.array([self.data[source]["ra"], self.data[source]["dec"]])
