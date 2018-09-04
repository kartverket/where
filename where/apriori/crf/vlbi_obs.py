"""A CRF Factory for creating sources with positions from VLBI observation files

Description:
------------

Reads source positions from the VLBI observation files in either NGS or vgosdb format.



"""

# External library imports
import numpy as np
import netCDF4

# Where imports
from where import apriori
from where.apriori import crf
from where.lib import config
from where.lib import log
from where.lib import plugins
from where import parsers
from where.lib import files
from where.lib import util


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
            return getattr(self, "_read_data_{}".format(fmt))()
        except AttributeError:
            log.fatal("Format '{}' is unknown for reference frame {}", fmt, self)

    def _read_data_vgosdb(self):
        """Read data from vgosdb observation files

        Returns:
            Dict:  Dictionary containing data about each source defined in this celestial reference frame.
        """
        # TODO make parser
        sources_file = files.path(file_key="vlbi_obs_vgosdb") / "Apriori" / "Source.nc"
        source_names = apriori.get("vlbi_source_names")
        d = netCDF4.Dataset(sources_file)
        sources = d["AprioriSourceList"][:]
        radec = d["AprioriSource2000RaDec"][:]
        return {
            source_names[str(s, "utf-8").strip()]["iers_name"]: dict(ra=coord[0], dec=coord[1])
            for s, coord in zip(sources, radec)
        }

    def _read_data_ngs(self):
        """Read data from NGS observation files

        Returns:
            Dict:  Dictionary containing data about each source defined in this celestial reference frame.
        """
        data = parsers.parse_key("vlbi_obs_ngs", parser="vlbi_ngs_sources").as_dict()
        src_names = apriori.get("vlbi_source_names")

        # Replace IVS name of source with official IERS name
        return {src_names[ivsname]["iers_name"]: coords for ivsname, coords in data.items() if ivsname in src_names}

    def _calculate_pos_crs(self, source):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            source (String):  Key specifying which source to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 2-vector.
        """
        return np.array([self.data[source]["ra"], self.data[source]["dec"]])
