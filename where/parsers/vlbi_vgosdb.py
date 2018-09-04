"""A parser for reading VLBI data from vgosDb files

Description:
------------

Reads data from files in the vgosDb files as defined in [1]. The data is organized in multiple smaller database
files based on netCDF.

References:
-----------

..[1] vgosDb format
    ftp://gemini.gsfc.nasa.gov/pub/misc/jmg/VLBI_Structure_2013Jun11.pdf

"""
# Standard library imports
import os

# External library imports
import netCDF4
import numpy as np
from scipy import interpolate

# Where imports
from where import apriori
from where.lib import config
from where.lib import constant
from where.lib import files
from where.lib import log
from where.lib.time import Time
from where.parsers import parser
from where.lib import plugins
from where.ext import sofa


@plugins.register
class VlbiVgosdbParser(parser.Parser):

    def __init__(self, rundate, session=None):
        super().__init__(rundate)
        self.file_key = "vlbi_obs_vgosdb"

        # Read arc-length from config file
        self.arc_length = config.tech.arc_length.int

        # Set session variables
        self.vars["session"] = session
        self.session_dir = files.path(self.file_key, file_vars=self.vars)
        self.data_available = self.session_dir.exists()

    def read_data(self):
        self.read_apriori()
        self.read_observables()
        self.read_obs_derived()
        self.read_obs_edit()
        self.read_crossreference()
        self.read_stations()

    def read_apriori(self):
        # Read station apriori
        nc = self._read_nc_file("Apriori", "Station.nc")
        self.data["sta_name"] = np.array([str(sta, "utf-8").strip() for sta in nc["AprioriStationList"]])
        self.data["sta_xyz"] = nc["AprioriStationXYZ"][:]
        log.info("Found stations: {}", ", ".join(self.data["sta_name"]))

        # Read source apriori

        source_names = apriori.get("vlbi_source_names")
        nc = self._read_nc_file("Apriori", "Source.nc")

        # import IPython; IPython.embed()
        self.data["source_names"] = np.array(
            [source_names[str(src, "utf-8").strip()]["iers_name"] for src in nc["AprioriSourceList"]]
        )
        self.data["source_ra"] = nc["AprioriSource2000RaDec"][:, 0]
        self.data["source_dec"] = nc["AprioriSource2000RaDec"][:, 1]

    def read_crossreference(self):
        nc = self._read_nc_file("CrossReference", "ObsCrossRef.nc")
        self.data["obs2baseline"] = self._get_nc_value(nc, "Obs2Baseline")
        self.data["obs2scan"] = self._get_nc_value(nc, "Obs2Scan")

        nc = self._read_nc_file("CrossReference", "SourceCrossRef.nc")
        self.data["scan2src"] = self._get_nc_value(nc, "Scan2Source")

    def read_observables(self):
        nc = self._read_nc_file("Observables", "TimeUTC.nc")
        self.data["time"] = [
            "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{}".format(*ymdhm, s)
            for ymdhm, s in zip(nc["YMDHM"][:], nc["Second"][:])
        ]
        self.data["num_obs"] = len(self.data["time"])

        nc = self._read_nc_file("Observables", "GroupDelay_bX.nc")
        self.data.setdefault("obs", {}).update(observed_delay_ferr=nc["GroupDelaySig"][:] * constant.c)

    def read_obs_derived(self):
        nc = self._read_nc_file("ObsDerived", "Cal-SlantPathIonoGroup_bX.nc")
        self.data.setdefault("obs", {}).update(iono_delay=nc["Cal-SlantPathIonoGroup"][:, 0] * constant.c)
        self.data.setdefault("obs", {}).update(iono_delay_ferr=nc["Cal-SlantPathIonoGroupSigma"][:, 0] * constant.c)
        try:
            self.data.setdefault("obs", {}).update(iono_quality=nc["Cal-SlantPathIonoGroupDataFlag"][:])
        except IndexError:
            # The iono flag is sometimes missing from the data
            self.data.setdefault("obs", {}).update(iono_quality=np.zeros(self.data["num_obs"]))
            log.warn("Missing quality flag for ionosphere")

    def read_obs_edit(self):
        nc = self._read_nc_file("ObsEdit", "GroupDelayFull_iIVS_bX.nc")
        #        nc = self._read_nc_file('ObsEdit', 'GroupDelayFull_bX.nc')
        self.data.setdefault("obs", {}).update(observed_delay=nc["GroupDelayFull"][:] * constant.c)

        try:
            nc = self._read_nc_file("ObsEdit", "Edit_iIVS_V004.nc")
            self.data.setdefault("obs", {}).update(data_quality=self._get_nc_value(nc, "DelayFlag"))
        except OSError:
            # The file is sometimes missing
            self.data.setdefault("obs", {}).update(data_quality=np.zeros(self.data["num_obs"]))
            log.warn("Missing data quality information")

    def read_stations(self):
        for sta in self.data["sta_name"]:
            nc = self._read_nc_file(sta, "TimeUTC.nc")
            sta_time = [
                "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{}".format(*ymdhm, s)
                for ymdhm, s in zip(nc["YMDHM"][:], nc["Second"][:])
            ]
            self.data.setdefault(sta, {}).update(time=sta_time)
            self.data.setdefault(sta, {}).update(num_obs=len(sta_time))

            if len(sta_time) < 5:
                log.warn("Only {} observations for {}", len(sta_time), sta)

            try:
                nc = self._read_nc_file(sta, "Cal-Cable.nc")
                self.data.setdefault(sta, {}).update(cable=nc["Cal-Cable"][:])
            except OSError:
                # Cable calibration is missing
                self.data.setdefault(sta, {}).update(cable=np.zeros(len(sta_time)))
                log.warn("Missing cable calibration information for {}", sta)

            try:
                nc = self._read_nc_file(sta, "Met.nc")
                self.data.setdefault(sta, {}).update(temperature=nc["TempC"][:])
                self.data.setdefault(sta, {}).update(pressure=nc["AtmPres"][:])
            except OSError:
                # Met data is missing
                self.data.setdefault(sta, {}).update(temperature=np.full(len(sta_time), np.nan))
                self.data.setdefault(sta, {}).update(pressure=np.full(len(sta_time), np.nan))
                log.warn("Missing meterological data for {}", sta)

    def write_to_dataset(self, dset):
        """Store VLBI data in a dataset

        Args:
            dset: The Dataset where data are stored.
        """
        dset.num_obs = len(self.data["time"])

        # Vgosdb indicies
        scan2src = self.data["scan2src"] - 1
        obs2scan = self.data["obs2scan"] - 1
        obs2baseline = self.data["obs2baseline"] - 1

        # Time
        dset.add_time("time", val=self.data["time"], scale="utc", format="isot")

        # Observations
        for obs_type, values in self.data["obs"].items():
            dset.add_float(obs_type, val=values)

        # Sources
        icrf = apriori.get("crf", session=dset.dataset_name)
        for i, src in enumerate(self.data["source_names"]):
            # Read source coordinates from ICRF file.
            if src in icrf:
                self.data["source_ra"][i] = icrf[src].pos.crs[0]
                self.data["source_dec"][i] = icrf[src].pos.crs[1]
                log.debug(
                    "Using coordinates ({:0.3f}, {:0.3f}) for {}", icrf[src].pos.crs[0], icrf[src].pos.crs[1], src
                )
            else:
                log.info(
                    "{} not found in apriori crf. Using coordinates ({:0.3f}, {:0.3f}) from vgosdb",
                    src,
                    self.data["source_ra"][i],
                    self.data["source_dec"][i],
                )

        dset.add_text("source", val=self.data["source_names"][scan2src[obs2scan]])
        dset.add_direction(
            "src_dir", ra=self.data["source_ra"][scan2src[obs2scan]], dec=self.data["source_dec"][scan2src[obs2scan]]
        )

        # Stations
        dset.add_text("station_1", val=self.data["sta_name"][obs2baseline][:, 0])
        dset.add_text("station_2", val=self.data["sta_name"][obs2baseline][:, 1])
        baseline = np.core.defchararray.add(np.core.defchararray.add(dset.station_1, "/"), dset.station_2)
        dset.add_text("baseline", val=baseline)
        dset.add_text("pass", val=np.core.defchararray.add(np.core.defchararray.add(baseline, "/"), dset.source))

        trf = apriori.get("trf", time=dset.time)
        station_codes = apriori.get("vlbi_station_codes")

        for site, xyz in zip(self.data["sta_name"], self.data["sta_xyz"]):
            if site in station_codes:
                cdp = station_codes[site]["cdp"]
                trf_site = trf[cdp]
            else:
                trf_site = trf.closest(xyz, max_distance=5)
                cdp = trf_site.key

            self.data["pos_" + site] = trf_site.pos.itrs
            log.debug("Using position {} for {}", np.mean(self.data["pos_" + site], axis=0), site)

            ivsname = station_codes[cdp]["ivsname"]
            domes = station_codes[cdp]["domes"]
            marker = station_codes[cdp]["marker"]
            self.data["station_" + site] = dict(cdp=cdp, domes=domes, ivsname=ivsname, marker=marker, site_id=cdp)

        # Positions
        itrs_pos_1 = np.array([self.data["pos_" + s][i, :] for i, s in enumerate(dset.station_1)])
        itrs_vel_1 = np.zeros((dset.num_obs, 3))

        dset.add_posvel(
            "site_pos_1",
            time="time",
            other="src_dir",
            itrs=np.concatenate((itrs_pos_1, itrs_vel_1), axis=1),
            write_level="operational",
        )
        itrs_pos_2 = np.array([self.data["pos_" + s][i, :] for i, s in enumerate(dset.station_2)])
        itrs_vel_2 = np.zeros((dset.num_obs, 3))
        dset.add_posvel(
            "site_pos_2",
            time="time",
            other="src_dir",
            itrs=np.concatenate((itrs_pos_2, itrs_vel_2), axis=1),
            write_level="operational",
        )

        # Station data
        sta_fields = set().union(*[v.keys() for k, v in self.data.items() if k.startswith("station_")])
        for field in sta_fields:
            dset.add_text(field + "_1", val=[self.data["station_" + s][field] for s in dset.station_1])
            dset.add_text(field + "_2", val=[self.data["station_" + s][field] for s in dset.station_2])

        dset.add_float("cable_delay_1", val=self._create_station_array(dset, "cable", 1) * constant.c)
        dset.add_float("cable_delay_2", val=self._create_station_array(dset, "cable", 2) * constant.c)
        dset.add_float("pressure_1", val=self._create_station_array(dset, "pressure", 1))
        dset.add_float("pressure_2", val=self._create_station_array(dset, "pressure", 2))
        dset.add_float("temperature_1", val=self._create_station_array(dset, "temperature", 1))
        dset.add_float("temperature_2", val=self._create_station_array(dset, "temperature", 2))

        # Station meta
        station_keys = sorted([k for k, v in self.data.items() if k.startswith("station_")])
        pos_keys = sorted([k for k, v in self.data.items() if k.startswith("pos_")])

        for sta_key, pos_key in zip(station_keys, pos_keys):
            sta_name = sta_key.split("_")[-1]
            cdp = self.data[sta_key]["cdp"]
            ivsname = station_codes[cdp]["ivsname"]
            longitude, latitude, height, _ = sofa.iau_gc2gd(2, self.data[pos_key][0, :])
            dset.add_to_meta(ivsname, "cdp", cdp)
            dset.add_to_meta(ivsname, "site_id", cdp)
            dset.add_to_meta(ivsname, "domes", station_codes[cdp]["domes"])
            dset.add_to_meta(ivsname, "marker", station_codes[cdp]["marker"])
            dset.add_to_meta(ivsname, "description", station_codes[cdp]["description"])
            dset.add_to_meta(ivsname, "longitude", longitude)
            dset.add_to_meta(ivsname, "latitude", latitude)
            dset.add_to_meta(ivsname, "height", height)
            if sta_name != ivsname:
                dset.add_to_meta(sta_name, "cdp", cdp)
                dset.add_to_meta(sta_name, "site_id", cdp)
                dset.add_to_meta(sta_name, "domes", station_codes[cdp]["domes"])
                dset.add_to_meta(sta_name, "marker", station_codes[cdp]["marker"])
                dset.add_to_meta(sta_name, "description", station_codes[cdp]["description"])
                dset.add_to_meta(sta_name, "longitude", longitude)
                dset.add_to_meta(sta_name, "latitude", latitude)
                dset.add_to_meta(sta_name, "height", height)

        dset.meta["tech"] = "vlbi"
        dset.add_to_meta("input", "type", "VGOSDB")
        dset.add_to_meta("input", "file", files.path(self.file_key, file_vars=self.vars).stem)

        master = apriori.get("vlbi_master_file")
        dset.meta["master_file"] = master.data.get((self.rundate.timetuple().tm_yday, self.vars["session"]), {})

    def _read_nc_file(self, directory, file_name):
        full_path = os.path.join(self.session_dir, directory, file_name)
        self.dependencies.append(full_path)
        log.debug("Read netCDF-data from {}", full_path)

        return netCDF4.Dataset(full_path)

    def _get_nc_value(self, nc, name):
        var = nc[name]
        if hasattr(var, "REPEAT"):
            if len(var[:]) > 1:
                return np.tile(var[:], var.REPEAT).reshape(var.REPEAT, -1)
            else:
                return np.tile(var[:], var.REPEAT)
        else:
            return nc[name][:]

    def _create_station_array(self, dset, field, suffix):
        dset.default_field_suffix = suffix
        array = np.zeros(dset.num_obs)

        for sta in self.data["sta_name"]:
            sta_idx = dset.filter(station=sta)
            sta_time = Time(self.data[sta]["time"]).utc.mjd
            func = interpolate.interp1d(
                sta_time,
                self.data[sta][field],
                bounds_error=False,
                fill_value=(self.data[sta][field][0], self.data[sta][field][-1]),
                assume_sorted=True,
            )
            array[sta_idx] = func(dset.time.utc.mjd[sta_idx])

        return array
