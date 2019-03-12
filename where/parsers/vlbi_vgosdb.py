"""A parser for reading VLBI data from vgosDb files

Description:
------------

Reads data from files in the vgosDb files as defined in [1]. The data is organized in multiple smaller database
files based on netCDF.

References:
-----------

..[1] vgosDb format
    ftp://gemini.gsfc.nasa.gov/pub/misc/jmg/VLBI_Structure_2013Jun11.pdf

    new url needed

"""
# Standard library imports
from datetime import datetime

# External library imports
import numpy as np
from scipy import interpolate

# Midgard imports
from midgard.dev import plugins

# Where imports
from midgard.math.constant import constant
from where.lib import log
from where.lib import files
from where.parsers._parser import Parser
from where import parsers


@plugins.register
class VgosDbParser(Parser):
    """A parser for reading VLBI data from a VGOS database
    """

    _STATION_FIELDS = {
        "temperature": {"filestub": "Met", "variable": "TempC", "factor": 1, "nan_value": -999},
        "pressure": {"filestub": "Met", "variable": "AtmPres", "factor": 1, "nan_value": -999},
        "cable_delay": {"filestub": "Cal-Cable", "variable": "Cal-Cable", "factor": constant.c, "nan_value": np.nan},
    }

    def __init__(self, file_path, encoding=None):
        super().__init__(file_path, encoding=None)
        self.raw = {}

    def read_data(self):
        """Parse the vgosdb wrapper file

        self.data will be populated with information from the netcdf files
        """
        with files.open_path(self.file_path, mode="rt") as fid:
            self._parse_file(fid)

        self._organize_data()

    def _parse_file(self, fid):
        for line in fid:
            if not line or line.startswith("!"):
                continue
            line = line.split()
            if "begin" in line[0].lower():
                self._parse_block(fid, line[1], name=" ".join(line[2:]))

    def _parse_block(self, fid, block, name="", directory=""):
        # print("Parsing {} {}".format(block, name))
        for line in fid:
            if not line or line.startswith("!"):
                continue
            line = line.split()
            if line[0].lower().startswith("end") and line[1] == block:
                # print("Finished {} {}".format(block, name))
                return
            elif line[0].lower().startswith("begin"):
                # recursive call
                self._parse_block(fid, line[1], name=" ".join(line[2:]), directory=directory)
            elif line[0].lower().startswith("default_dir"):
                directory = line[1]
            elif line[0].endswith(".nc"):
                file_path = self.file_path.parents[0] / directory / line[0]
                if directory:
                    data = self.raw.setdefault(directory, {})
                else:
                    data = self.raw
                nc_name = file_path.stem.split("_")
                nc_stub = nc_name.pop(0)
                data = data.setdefault(nc_stub, {})
                for part in nc_name:
                    if part.startswith("b"):
                        data = data.setdefault(part[1:], {})

                # print("Parse {}".format(file_path))
                netcdf_data = parsers.parse_file("vlbi_netcdf", file_path=file_path).as_dict()
                if "TimeUTC" in file_path.stem:
                    self._parse_time(netcdf_data)
                data.update(netcdf_data)
            else:
                data = self.raw.setdefault(block, {})
                if name:
                    data = data.setdefault(name, {})
                data[line[0]] = " ".join(line[1:])

    def _organize_data(self):
        """ Copy content from self.raw to self.data and convert all data to arrays with num_obs length
        """
        meta = self.data.setdefault("meta", {})
        meta["session_code"] = self.raw["Session"].get("Session")

        # Epoch info
        self.data["time"] = self.raw["Observables"]["TimeUTC"]["time"]

        num_obs = len(self.data["time"])
        self.data["station_1"] = self.raw["Observables"]["Baseline"]["Baseline"].reshape(num_obs, -1)[:, 0]
        self.data["station_2"] = self.raw["Observables"]["Baseline"]["Baseline"].reshape(num_obs, -1)[:, 1]
        self.data["source"] = self.raw["Observables"]["Source"]["Source"]

        # Obs info
        try:
            self.data["observed_delay_ferr"] = self.raw["Observables"]["GroupDelay"]["X"]["GroupDelaySig"] * constant.c
        except KeyError:
            self.data["observed_delay_ferr"] = np.zeros(num_obs)
            log.error("Missing group delay formal error information")

        try:
            self.data["data_quality"] = self.raw["ObsEdit"]["Edit"]["DelayFlag"]
        except KeyError:
            self.data["data_quality"] = np.full(num_obs, np.nan)
            log.warn("Missing data quality information")

        try:
            self.data["observed_delay"] = self.raw["ObsEdit"]["GroupDelayFull"]["X"]["GroupDelayFull"] * constant.c
        except KeyError:
            self.data["observed_delay"] = np.full(num_obs, np.nan)
            log.error("Missing full group delay information")

        try:
            self.data["iono_delay"] = (
                self.raw["ObsDerived"]["Cal-SlantPathIonoGroup"]["X"]["Cal-SlantPathIonoGroup"].reshape(num_obs, -1)[
                    :, 0
                ]
                * constant.c
            )
        except KeyError:
            self.data["iono_delay"] = np.full(num_obs, np.nan)
            log.warn("Missing ionosphere delay information")

        try:
            self.data["iono_delay_ferr"] = (
                self.raw["ObsDerived"]["Cal-SlantPathIonoGroup"]["X"]["Cal-SlantPathIonoGroupSigma"].reshape(
                    num_obs, -1
                )[:, 0]
                * constant.c
            )
        except KeyError:
            self.data["iono_delay_ferr"] = np.zeros(len(self.data["time"]))
            log.warn("Missing ionosphere delay formal error information")

        try:
            self.data["iono_quality"] = self.raw["ObsDerived"]["Cal-SlantPathIonoGroup"]["X"][
                "Cal-SlantPathIonoGroupDataFlag"
            ]
        except KeyError:
            self.data["iono_quality"] = np.full(len(self.data["time"]), np.nan)
            log.warn("Missing ionosphere quality information")

        # Station dependent info
        for field, params in self._STATION_FIELDS.items():
            self.data[field + "_1"] = np.zeros(len(self.data["time"]))
            self.data[field + "_2"] = np.zeros(len(self.data["time"]))
            for station in self.raw["Head"]["StationList"]:
                sta_idx_1 = self.data["station_1"] == station
                sta_idx_2 = self.data["station_2"] == station
                sta_key = station.replace(" ", "_")
                sta_time = self.raw[sta_key]["TimeUTC"]["sec_since_ref"]
                try:
                    sta_data = self.raw[sta_key][params["filestub"]][params["variable"]]
                    missing_idx = np.isclose(sta_data, params["nan_value"])
                    sta_data[missing_idx] = np.nan
                    if missing_idx.any():
                        log.warn(f"Missing {field} data for {station}")
                except KeyError:
                    sta_data = np.full(len(sta_time), np.nan)
                    log.warn(f"Missing all {field} data for {station}")

                if len(sta_data) == 1:
                    # Use constant function if there is only one data point
                    func = lambda _: sta_data[0]
                else:
                    func = interpolate.interp1d(
                        sta_time,
                        sta_data,
                        bounds_error=False,
                        fill_value=(sta_data[0], sta_data[-1]),
                        assume_sorted=True,
                    )
                epochs_1 = self.raw["Observables"]["TimeUTC"]["sec_since_ref"][sta_idx_1]
                epochs_2 = self.raw["Observables"]["TimeUTC"]["sec_since_ref"][sta_idx_2]
                self.data[field + "_1"][sta_idx_1] = func(epochs_1) * params["factor"]
                self.data[field + "_2"][sta_idx_2] = func(epochs_2) * params["factor"]

    def _parse_time(self, time_dict):
        part1 = time_dict.pop("YMDHM")
        part2 = time_dict.pop("Second")
        time_dict["time"] = ["{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{}".format(*t1, t2) for t1, t2 in zip(part1, part2)]

        self.raw["dt_0"] = datetime(*part1[0], int(part2[0]))
        time_dict["sec_since_ref"] = np.array(
            [(datetime(*t1, int(t2)) - self.raw["dt_0"]).total_seconds() for t1, t2 in zip(part1, part2)]
        )
