"""A TRF Factory for creating sites with positions from the NMA database

Description:
------------

Reads station positions through a web service. A time series of positions is available, so positions are given
according to the time series.




"""
# Standard library imports
from datetime import datetime, timedelta
import json

# External library imports
import numpy as np
import requests

# Where imports
from where.apriori import trf
from where.lib import config
from where.lib import log
from where.lib import plugins


@plugins.register
class NmaTrf(trf.TrfFactory):
    """A factory for using positions from the custom TRF-config file.
    """

    @property
    def url(self):
        """URLs used to access the database
        """
        return dict(site=config.where.database.nma_ws_site.str, sites=config.where.database.nma_ws_sites.str)

    @property
    def sites(self):
        url = self.url["sites"]
        sites = json.loads(requests.get(url).text)
        return [c["fourCharId"] for s in sites for c in s["siteConfigs"]]

    def site(self, key):
        """Positions and information about one site in the reference frame

        Args:
            key (String):  Key specifying which site to calculate position for.

        Returns:
            TrfSite:  Object with positions and information about site.
        """
        if key not in self.data:
            url = self.url["site"].format(site=key)
            log.info("Reading information about {} from {}", key, url)
            db_data = json.loads(requests.get(url).text)
            if not db_data:
                log.warn("No information returned for {}", key)
                self._data[key] = dict(id=None, provider=None, siteConfig_id=None)
                return super().site(key)

            site_fields = ("id", "provider", "siteConfig_id")
            site_data = {f: db_data[-1][f] for f in site_fields}

            epoch_data = list()
            epoch_fields = ("year", "doyStart", "doyEnd", "x", "y", "z")
            for epoch in db_data:
                if self.version is not None and epoch["geodeticDatum"]["geodeticDatumName"] != self.version:
                    continue
                epoch_data.append([epoch[f] for f in epoch_fields])
            epoch_data = np.array(epoch_data)

            site_data["time_start"] = np.array(
                [datetime.strptime("{:02.0f} {:.0f}".format(y, d), "%y %j") for y, d in epoch_data[:, [0, 1]]]
            )
            site_data["time_end"] = np.array(
                [datetime.strptime("{:02.0f} {:.0f}".format(y, d), "%y %j") for y, d in epoch_data[:, [0, 2]]]
            ) + timedelta(
                days=1
            )
            site_data["pos"] = epoch_data[:, [3, 4, 5]]
            self._data[key] = site_data

        return super().site(key)

    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Data are stored in the database, instead of reading it all we provide an empty dictionary that will be filled
        as sites are requested.

        Returns:
            Dict:  Empty dictionary that will contain data about each site requested.
        """
        return dict()

    def _calculate_pos_itrs(self, site):
        """Calculate positions for the given time epochs

        There are no velocities available, so same position is returned for all time epochs

        Args:
            site (String):  Key specifying which site to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        site_data = self.data[site]
        pos = np.full((len(self.time), 3), np.nan)
        if "pos" not in site_data:
            return pos

        for idx, epoch in enumerate(self.time.datetime):
            try:
                ep_idx = np.where((site_data["time_start"] <= epoch) & (epoch < site_data["time_end"]))[0][0]
            except IndexError:
                continue
            pos[idx] = site_data["pos"][ep_idx]
        return pos
