"""A TRF Factory for creating sites with positions from a custom config file

Description:
------------

Reads station positions and velocities from a custom config file. No velocities are available in the file, so the same
position is used for all time epochs.



$Revision: 15220 $
$Date: 2018-05-31 13:51:15 +0200 (Thu, 31 May 2018) $
$LastChangedBy: hjegei $
"""
# External library imports
import numpy as np

# Where imports
from midgard.config.config import Configuration
from where.apriori import trf
from where.lib import config
from where.lib import dependencies
from where.lib import files
from where.lib import plugins


@plugins.register
class CustomTrf(trf.TrfFactory):
    """A factory for using positions from the custom TRF-config file.
    """

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Delegates to _read_data_<self.version> to read the actual data.

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        trf = Configuration("trf_custom")
        trf.profiles = config.analysis.get("tech", value=self.version, default="").list
        trf_path = files.path("trf-custom")
        trf_local_path = files.path("trf-custom_local")

        trf.update_from_file(trf_path)
        dependencies.add(trf_path)
        if trf_local_path.exists():
            trf.update_from_file(trf_local_path)
            dependencies.add(trf_local_path)

        data = dict()
        for key in trf.sections:
            info = {k: v for k, v in trf[key].as_dict().items() if not k == "pos_itrs"}
            info["pos"] = np.array(trf[key].pos_itrs.list, dtype=float)
            data[key] = info

        return data

    def _calculate_pos_itrs(self, site):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            site (String):  Key specifying which site to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        if self.time.size == 1:
            return self.data[site]["pos"]
        return self.data[site]["pos"][None, :].repeat(self.time.size, axis=0)
