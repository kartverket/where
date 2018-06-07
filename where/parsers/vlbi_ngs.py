"""A parser for reading VLBI data from NGS files

Description:
------------

Reads data from files in the NGS file format as defined in http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt
(revision date June 11, 2007) [1].

References:
-----------

[1] NGS file format.
    http://lacerta.gsfc.nasa.gov/mk5/help/dbngs_format.txt

[2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010).
    http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html



"""

# Standard library imports
import itertools
import os.path
from datetime import timedelta

# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import config
from where.lib import constant
from where.lib import files
from where.lib import log
from where.parsers import parser
from where.lib import plugins
from where.ext import sofa
from where.lib.unit import unit


@plugins.register
class VlbiNgsParser(parser.Parser):
    """A parser for reading VLBI data from NGS files
    """

    def __init__(self, rundate, session):
        super().__init__(rundate)
        self.file_key = "vlbi_obs_ngs"
        self.vars["session"] = session

        # Read arc-length from config file
        self.arc_length = config.tech.arc_length.int

    def read_data(self):
        """Read input data from all session files within the arc length
        """
        date_to_read = self.rundate
        while date_to_read < self.rundate + timedelta(days=self.arc_length):
            self.vars.update(config.date_vars(date_to_read))
            date_to_read += timedelta(days=1)

            # Find latest version if necessary
            if not os.path.exists(files.path(self.file_key, file_vars=self.vars)):
                versions = files.glob_variable(self.file_key, "obs_version", r"\d{3}", file_vars=self.vars)
                try:
                    self.vars["obs_version"] = sorted(versions)[-1]
                except IndexError:
                    self.vars["obs_version"] = "xxx"  # Used for version in filename in warning about missing data
                    continue

            # Read and parse file
            self.dependencies.append(files.path(self.file_key, file_vars=self.vars))
            with files.open(self.file_key, file_vars=self.vars, mode="rt") as fid:
                self.parse_file(fid)

        if not self.data:
            self.data_available = False

    #
    # PARSERS for reading each line of the NGS file.
    #
    def setup_parsers(self):
        # Ignore the header, first two lines
        header_parser = parser.define_parser(
            end_marker=lambda _l, line_num, _n: line_num == 2, label=None, parser_def=None
        )

        # Each line defines a station (site)
        station_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line == "$END",
            label=lambda line, _ln: line != "$END" and "station",
            parser_def={
                "station": {
                    "parser": self.parse_station,
                    "fields": {"name": (0, 8), "pos_x": (10, 25), "pos_y": (25, 40), "pos_z": (40, 55)},
                }
            },
        )

        # Each line defines a radio source
        source_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line == "$END",
            label=lambda line, _ln: line != "$END" and "source",
            parser_def={
                "source": {
                    "parser": self.parse_radio_source,
                    "fields": {
                        "name": (0, 8),
                        "ra_hrs": (10, 12),
                        "ra_mins": (13, 15),
                        "ra_secs": (16, 28),
                        "dec_degs": (29, 32),
                        "dec_mins": (33, 35),
                        "dec_secs": (36, 48),
                    },
                }
            },
        )

        # One line with auxiliary parameters
        param_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line == "$END",
            label=lambda line, _ln: line != "$END" and "param",
            parser_def={
                "param": {
                    "parser": self.parse_session,
                    "fields": {
                        "ref_freq": (0, 20),
                        "group_delay": (20, 30),
                        "delay_type": (31, 33),
                        "delay_rate_type": (34, 36),
                    },
                }
            },
        )

        # Observations are listed on 9 lines
        obs_parser = parser.define_parser(
            end_marker=lambda _l, _ln, next_line: next_line[78:80] == "01",
            end_callback=self.copy_cache_to_obs,
            label=lambda line, _ln: line[78:80],
            parser_def={
                "01": {
                    "parser": self.parse_obs_meta,
                    "fields": {
                        "station_1": (0, 8),
                        "station_2": (10, 18),
                        "source": (20, 28),
                        "year": (29, 33),
                        "month": (34, 36),
                        "day": (37, 39),
                        "hour": (40, 42),
                        "minute": (43, 45),
                        "seconds": (46, 60),
                    },
                },
                "02": {
                    "parser": self.parse_obs(unit_in="nanoseconds", except_fields=("data_quality",)),
                    "fields": {
                        "observed_delay": (0, 20),
                        "observed_delay_ferr": (20, 30),
                        #                                  'observed_delay_rate':      (30, 50),
                        #                                  'observed_delay_rate_ferr': (50, 60),
                        "data_quality": (60, 62),
                        #                                  'flag_delay_type':          (63, 65),
                        #                                  'flag_delay_rate_type':     (66, 68),
                    },
                },
                "03": {
                    "parser": self.parse_obs(),
                    "fields": {
                        "correlation_coeff": (0, 10),
                        "correlation_coeff_ferr": (10, 20),
                        "fringe_amplitude": (20, 30),
                        "fringe_amplitude_ferr": (30, 40),
                        "total_fringe_phase": (40, 60),
                        "total_fringe_phase_ferr": (60, 70),
                    },
                },
                "05": {
                    "parser": self.parse_obs(unit_in="nanoseconds"),
                    "fields": {"cable_delay_1": (0, 10), "cable_delay_2": (10, 20)},
                },
                "06": {
                    "parser": self.parse_obs_missing,
                    "fields": {
                        "temperature_1": (0, 10),
                        "temperature_2": (10, 20),
                        "pressure_1": (20, 30),
                        "pressure_2": (30, 40),
                    },
                },
                "08": {
                    "parser": self.parse_obs(unit_in="nanoseconds", except_fields=("iono_quality",)),
                    "fields": {
                        "iono_delay": (0, 20),
                        "iono_delay_ferr": (20, 30),
                        #                                  'iono_delay_rate':           (30, 50),
                        #                                  'iono_delay_rate_ferr':      (50, 60),
                        "iono_quality": (61, 63),
                    },
                },
            },
        )

        return itertools.chain(
            [header_parser, station_parser, source_parser, param_parser], itertools.repeat(obs_parser)
        )

    def parse_station(self, line, _):
        """Read station position

        Reads the station position from the NGS file. There is no consistent way to link stations listed in the NGS
        file with stations in ITRF apriori data, so use the position and look up the closest VLBI station in ITRF
        within a radius of 5 meters.

        Stores both a 3-vector XYZ-position, and a 3-vector LatLongHeight-position. For the transformation to
        LatLongHeight, the IERS conventions recommend using the GRS80 ellipsoid [2, section 4.2.6].

        Args:
            line:  Input data from NGS file

        """
        self.data["pos_" + line["name"]] = np.array([float(line["pos_x"]), float(line["pos_y"]), float(line["pos_z"])])

    def parse_radio_source(self, line, _):
        """Read radio source coordinates

        Reads the radio source coordinates from the ICRF and NGS files. Stores right ascension and declination in
        radians.

        Args:
            line:  Input data from NGS file
        """
        ra = unit.hms_to_rad(float(line["ra_hrs"]), int(line["ra_mins"]), float(line["ra_secs"]))
        dec = unit.dms_to_rad(float(line["dec_degs"].replace(" ", "")), int(line["dec_mins"]), float(line["dec_secs"]))
        self.data["src_" + line["name"]] = np.array([ra, dec])

    def parse_session(self, line, _):
        """Read session information

        Args:
            line:   Input data from NGS file.
        """
        self.data.setdefault("session", {}).update(line)

    def parse_obs_meta(self, line, cache):
        """Reads meta information like time stamp and station and source id

        Creates a new observation on the dataset with proper metainformation. Observation data are added to the dataset
        later by a calculator as it is more effective not to continuously resize numpy arrays.

        Args:
            line:  Input data from NGS file

        """
        line["seconds"] = "{:013.10f}".format(float(line["seconds"]))
        try:
            line["year"] = "{:4d}".format(int(line["year"]))
        except ValueError as e:
            # In some sessions the year field is '19 0' when it should be '2000'
            year = line["year"].split()
            if len(year) == 2:
                if int(year[1]) < 50:
                    line["year"] = "{:4d}".format(2000 + int(year[1]))
                else:
                    line["year"] = "{:4d}".format(1900 + int(year[1]))
            else:
                raise e

        obs = {
            "time": "{year:0>4}-{month:0>2}-{day:0>2}T{hour:0>2}:{minute:0>2}:{seconds}".format(**line),
            "station_1": line["station_1"],
            "station_2": line["station_2"],
            "source": line["source"].replace(".", "dot"),
            "pass": "{station_1}/{station_2}/{source}".format(**line),
            "baseline": "{station_1}/{station_2}".format(**line),
        }
        cache["obs"] = obs

        for field, value in obs.items():
            self.meta.setdefault(field, list()).append(value)

    def parse_obs(self, unit_in="meter", except_fields=()):
        """Read information about an observation

        Stores the information in the temporary cache-dict, which will be transfered to self.data when all information
        about this observation is parsed. If `unit_in` is given, all values will be converted to meter unless the field
        is listed in the `except_fields`-list.

        Args:
            unit_in (String):      Name of unit of values to be parsed.
            except_fields (Tuple): Names of fields where values should not be converted to meters.
        """
        # Find scale factor for converting to meter
        if unit_in == "meter":
            scale_factor = 1
        else:
            quantity = unit(unit_in)
            try:
                scale_factor = quantity.to("meter").magnitude
            except unit.DimensionalityError:
                # Try to convert between time and length by multiplying by the speed of light
                scale_factor = (quantity * constant.c * unit("meters per second")).to("meter").magnitude

        # Define the function doing the actual parsing
        def parse_func(line, cache):
            for field in line:
                if field.startswith("flag_"):  # Flags are currently ignored
                    continue
                try:
                    cache.setdefault("values", {})[field] = float(line[field])
                except ValueError:
                    cache.setdefault("values", {})[field] = 0
                    log.debug(
                        "Could not convert {} = {} to a number for {} {}. Value set to 0.0.",
                        field,
                        line[field],
                        cache["obs"]["pass"],
                        cache["obs"]["time"],
                    )
                if field not in except_fields:
                    cache["values"][field] *= scale_factor

        return parse_func

    def parse_obs_missing(self, line, cache):
        for key, value in line.items():
            if value.startswith("-999"):
                line[key] = "nan"
        self.parse_obs()(line, cache)

    def copy_cache_to_obs(self, cache):
        """Copy temporary cache storage to data.obs list of observation data

        Args:
            cache: Temporary storage for information
        """
        self.data.setdefault("obs", []).append(cache["values"])

    #
    # CALCULATORS for calculating other necessary data
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        if self.data:
            return []
        else:
            self.data_available = False
            return list()

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Store VLBI data in a dataset

        Args:
            dset: The Dataset where data are stored.
        """
        dset.num_obs = len(self.data["obs"])
        dset.add_time("time", val=self.meta.pop("time"), scale="utc", format="isot", write_level="operational")
        for field, value in self.meta.items():
            dset.add_text(field, val=value, write_level="operational")

        # Observations
        obs_fields = set().union(*[d.keys() for d in self.data["obs"]])
        for field in obs_fields:
            dset.add_float(field, val=np.array([o[field] for o in self.data["obs"]]), write_level="operational")

        # Source directions
        icrf = apriori.get("crf", session=dset.dataset_name)
        ra = np.array([icrf[s].pos.crs[0] if s in icrf else 0 for s in dset.source])
        dec = np.array([icrf[s].pos.crs[1] if s in icrf else 0 for s in dset.source])

        # If there are more than 300 sources in a NGS-file the source names are gibberish
        bad_source_idx = np.where(ra == 0)[0]
        bad_sources = np.array(dset.source)[bad_source_idx]
        for s in np.unique(bad_sources):
            # TODO: discard automatically? obs_format dependent edit file is bad
            log.check("Unknown source {}. Add to vlbi_ignore_source", s)

        dset.add_direction("src_dir", ra=ra, dec=dec, write_level="operational")

        # Look up positions in ITRF file
        log.info("Found stations: {}", ", ".join(dset.unique("station")))
        trf = apriori.get("trf", time=dset.time)
        station_codes = apriori.get("vlbi_station_codes")
        for site in dset.unique("station"):
            if site in station_codes:
                cdp = station_codes[site]["cdp"]
                trf_site = trf[cdp]
            else:
                trf_site = trf.closest(self.data["pos_" + site], max_distance=5)
                cdp = trf_site.key
                log.warn("Undefined station name {}. Assuming station is {}.".format(site, trf_site.name))

            self.data["pos_" + site] = trf_site.pos.itrs
            log.debug(
                "Using position {} for {} from {}", np.mean(self.data["pos_" + site], axis=0), site, trf_site.source
            )

            ivsname = station_codes[cdp]["ivsname"]
            domes = station_codes[cdp]["domes"]
            marker = station_codes[cdp]["marker"]
            self.data["station_" + site] = dict(site_id=cdp, cdp=cdp, domes=domes, ivsname=ivsname, marker=marker)

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
            dset.add_text(
                field + "_1", val=[self.data["station_" + s][field] for s in dset.station_1]
            )  # write_level='analysis')
            dset.add_text(
                field + "_2", val=[self.data["station_" + s][field] for s in dset.station_2]
            )  # write_level='analysis')

        # Station meta
        station_keys = sorted([k for k, v in self.data.items() if k.startswith("station_")])
        pos_keys = sorted([k for k, v in self.data.items() if k.startswith("pos_")])

        for sta_key, pos_key in zip(station_keys, pos_keys):
            sta_name = sta_key.split("_")[-1]
            cdp = self.data[sta_key]["cdp"]
            ivsname = station_codes[cdp]["ivsname"]
            longitude, latitude, height, _ = sofa.iau_gc2gd(2, self.data[pos_key][0, :])  # TODO: Reference ellipsoid
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

        dset.add_to_meta("input", "file", files.path(self.file_key, file_vars=self.vars).stem)
        dset.add_to_meta("input", "type", "NGS")
        dset.meta["tech"] = "vlbi"

        master = apriori.get("vlbi_master_file")
        dset.meta.update(master.get((self.rundate.timetuple().tm_yday, self.vars["session"]), {}))
