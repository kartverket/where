"""A parser for reading SP3-d orbit files

Example:
--------

    from where import data
    from where import parsers

    # Parse data
    parser = parsers.parse('orbit_sp3d', file_path=file_path, rundate=rundate)

    # Create a empty Dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=0, empty=True)

    # Fill Dataset with parsed data
    parser.write_to_dataset(dset)


Description:
------------

Reads data from SP3-d files (see :cite:`hilla2016`). The reference frame is defined in the header of the SP3 file,
which is for IGS products the current IGS terrestrial reference frame (e.g. IGb08). The time system is for IGS products
the GPS time scale. The orbit position and velocities are given normally for every 15 minutes.

The file to be parsed should be specified like:

    from where import parsers
    parser = parsers.parse('orbit_sp3c', file_path=file_path, rundate=rundate)


TODO:
-----
Important: Bad satellites has to be removed from observation Dataset!!!! E.g. indicated by bad satellite clock values.



"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers import parser
from where.lib import config
from where.lib import constant
from where.lib import log
from where.lib.unit import unit


@plugins.register
class OrbitSp3dParser(parser.Parser):
    """A parser for reading data from SP3-d orbit files

    Attributes:
        data (dict):             Dict containing the (precise orbit) data read from file.
        data_available (bool):   Indicator of whether data are available.
        dependencies (list):     List of files that have been read by the parser.
        file_key (str):          Key to the SP3 orbit file defined in files.conf file.
        file_path (pathlib.PosixPath):  File path to SP3 orbit file.
        meta (dict):             Dict containing the metainformation read from file, whereby the current days define the
                                 keys and the values are the SP3 header entries.
        rundate (datetime.date): Date of model run.
        vars (dict):             Variables needed to define date identifiers for saving SP3 header information in 'meta'
                                 variable.


    Methods:
        calculate_data()        Carry out the defined calculators
        copy_cache_to_data()    Copy contents of the cache to the data datastructure
        copy_cache_to_meta()    Copy contents of the cache to the meta datastructure
        parse()                 Parse data
        parse_default()         Add the contents of line to data
        parse_default_meta()    Add the contents of line to meta
        parse_file()            Read a data file and parse the content
        parse_line()            Parse line
        read_data()             Read data from datafiles
        setup_calculators()     List steps necessary for postprocessing
        setup_parsers()         Setup parser definition
        write_to_dataset()      Write data based on GNSS SP3 orbit file

        _parse_date()           Parse date orbit position/velocity block
        _parse_float()          Parse float entries of SP3 header to instance variable 'meta'
        _parse_position()       Parse orbit position
        _parse_string()         Parse string entries of SP3 header to instance variable 'meta'
        _parse_velocity()       Parse orbit velocity
    """

    def __init__(self, rundate, file_path):
        """Set up the basic information needed by the parser

        Args:
            rundate (datetime.date):    The model run date.
            file_path (str):            Path to orbit-file to parse.
        """
        super().__init__(rundate=rundate, file_path=file_path)
        self.vars = config.date_vars(rundate)

    #
    # PARSER for reading each line of the SP3 file.
    #
    def setup_parsers(self):

        # Parser for SP3 header
        header_parser = parser.define_parser(
            end_marker=lambda _l, _ln, next_line: next_line.startswith("*"),
            # end_marker=lambda _l, line_num, _n: line_num == 22,
            label=lambda _l, line_num: line_num,
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6
                # #dP2016  3  1  0  0  0.00000000      96 ORBIT IGb08 HLM  IGS
                1: {
                    "parser": self._parse_string,
                    "fields": {
                        "version": (1, 2),
                        "pv_flag": (2, 3),
                        "year": (3, 7),
                        "month": (8, 10),
                        "day": (11, 13),
                        "hour": (14, 16),
                        "minute": (17, 19),
                        "second": (20, 31),
                        "num_epoch": (32, 39),
                        "data_used": (40, 45),
                        "coord_sys": (46, 51),
                        "orb_type": (52, 55),
                        "agency": (56, 60),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6
                # ## 1886 172800.00000000   900.00000000 57448 0.0000000000000
                2: {
                    "parser": self._parse_string,
                    "fields": {
                        "gpsweek": (3, 7),
                        "gpssec": (8, 23),
                        "epoch_interval": (24, 38),
                        "mjd_int": (39, 44),
                        "mjd_frac": (45, 60),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6
                # %c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
                15: {"parser": self._parse_string, "fields": {"file_type": (3, 5), "time_sys": (9, 12)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6
                # %f  1.2500000  1.025000000  0.00000000000  0.000000000000000
                17: {"parser": self._parse_float, "fields": {"base_posvel": (3, 13), "base_clkrate": (14, 26)}},
            },
        )

        # Parser for SP3 data block

        # *  2016  3  1  0  0  0.00000000
        # PG01  10138.887745 -20456.557725 -13455.830128     13.095853  7  6  4 137
        # PG02 -21691.921884  13338.131173  -6326.904893    599.417359 10  7  8 114
        # PG03   1061.483783 -15622.426751 -21452.532447    -30.693182  7  6 10 125
        # PG04  25398.213954   6966.881030   4188.487313 999999.999999
        # PG05 -20431.536257   5143.439387  16220.688245   -141.634044  9  8  7 101
        # PG06 -19896.425641    760.039334 -17564.616810    153.322562 12 10 13 104
        # PG07  -5499.233721 -14614.719047  21564.675152    467.390140  7 10  7 123
        # PG08   7359.468852 -20268.667422  15550.334189    -25.778533  8  6  8 113
        data_parser = parser.define_parser(
            end_marker=lambda _l, _ln, next_line: next_line.startswith("*"),
            label=lambda line, _ln: line[0],
            parser_def={
                "*": {
                    "parser": self._parse_date,
                    "fields": [None, "year", "month", "day", "hour", "minute", "second"],
                },
                "P": {
                    "parser": self._parse_position,
                    "fields": {
                        "sat": (1, 4),
                        "pos_x": (4, 18),
                        "pos_y": (18, 32),
                        "pos_z": (32, 46),
                        "clk_bias": (46, 60),
                        "sig_pos_x": (61, 63),
                        "sig_pos_y": (64, 66),
                        "sig_pos_z": (67, 69),
                        "sig_clk_bias": (70, 73),
                        "clk_event_flag": (74, 75),
                        "clk_pred_flag": (75, 76),
                        "maneuver_flag": (78, 79),
                        "orb_pred_flag": (79, 80),
                    },
                },
                "V": {
                    "parser": self._parse_velocity,
                    "fields": {
                        "sat": (1, 4),
                        "vel_x": (4, 18),
                        "vel_y": (18, 32),
                        "vel_z": (32, 46),
                        "clk_rate": (46, 60),
                        "sig_vel_x": (61, 63),
                        "sig_vel_y": (64, 66),
                        "sig_vel_z": (67, 69),
                        "sig_clk_rate": (70, 73),
                    },
                },
            },
        )

        return itertools.chain([header_parser], itertools.repeat(data_parser))

    def _parse_date(self, line, cache):
        """Parse date orbit position/velocity block

        Args:
            line (dict):  Dict containing the fields of a line.
            cache (dict): Temporary dictionary with the fields 'key' and 'values'.
        """
        # Create Time object
        cache["time"] = "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}" "".format(
            year=int(line["year"]),
            month=int(line["month"]),
            day=int(line["day"]),
            hour=int(line["hour"]),
            minute=int(line["minute"]),
            second=float(line["second"]),
        )

    def _parse_float(self, line, _):
        """Parse float entries of SP3 header to instance variable 'meta'

        Args:
            line (dict):  Dict containing the fields of a line.
        """
        date = "{year}-{month:02d}-{day:02d}".format(
            year=int(self.vars["yyyy"]), month=int(self.vars["m"]), day=int(self.vars["d"])
        )  # TODO: What would be a better solution for getting current date?

        for k, v in line.items():
            self.meta.setdefault(date, dict()).update({k: float(v)})

    def _parse_position(self, line, cache):
        """Parse orbit position

        Bad or absent position values indicated by 0.000000 and clock values indicated by 999999.999999 are set to
        not a number 'nan'.

        Satellite identifier (e.g. G01) consists of satellite system identicator 'G' and a satellite number '01'. In the
        SP3-a format the satellite system identicator is not used. SP3-a was designed only for GPS satellites, therefore
        if the format version 'a' is given, the satellite system identicator 'G' is introduced before the satellite
        number.

        The standard deviations of the satellite and the clock correction has to multiplied with the base floating point
        numbers defined in header line 15.

        Args:
            line (dict):  Dict containing the fields of a line.
            cache (dict): Temporary dictionary with the fields 'key' and 'values'.
        """

        # Remove identical epochs given in two different SP3 files
        if cache["line_num"] == 2:
            if "time" in self.data:
                if cache["time"] in self.data["time"]:
                    log.warn("Identical epoch {} given in the SP3 files.", cache["time"])
                    return

        date = "{year}-{month:02d}-{day:02d}".format(
            year=int(self.vars["yyyy"]), month=int(self.vars["m"]), day=int(self.vars["d"])
        )  # TODO: What would be a better solution for getting current date?

        # SP3-a (GPS-only) format files missing satellite system identicator 'G' before satellite number
        if self.meta[date]["version"] == "a":
            line["sat"] = "G" + line["sat"].zfill(2)

        # Set bad or absent positional values to 'nan' indicated by 0.000000
        for k in list(["pos_x", "pos_y", "pos_z"]):
            if float(line[k]) == 0.0:
                line[k] = float("nan")

        # Set bad or absent clock bias values to 'nan' indicated by 999999.999999
        if float(line["clk_bias"]) == 999999.999999:
            line["clk_bias"] = float("nan")

        # Set not given sigmas in SP3 file to 'nan'
        for k in list(["sig_pos_x", "sig_pos_y", "sig_pos_z", "sig_clk_bias"]):
            if line[k] == "":
                line[k] = float("nan")

        self.data.setdefault("time", list()).append(cache["time"])
        self.data.setdefault("satellite", list()).append(line["sat"])
        self.data.setdefault("sat_pos", list()).append(
            np.array([float(line["pos_x"]), float(line["pos_y"]), float(line["pos_z"])]) * unit.kilometer2meter
        )
        self.data.setdefault("sat_clock_bias", list()).append(
            float(line["clk_bias"]) * unit.microsecond2second * constant.c
        )
        self.data.setdefault("sat_pos_sigma", list()).append(
            np.array([float(line["sig_pos_x"]), float(line["sig_pos_y"]), float(line["sig_pos_z"])])
            * self.meta[date]["base_posvel"]
            * unit.millimeter2meter
        )
        self.data.setdefault("sat_clock_bias_sigma", list()).append(
            float(line["sig_clk_bias"]) * self.meta[date]["base_clkrate"] * unit.picosecond2second * constant.c
        )

        # Get GNSS identifier
        sys = line["sat"][0]
        self.data.setdefault("system", list()).append(sys)

    def _parse_string(self, line, _):
        """Parse string entries of SP3 header to instance variable 'meta'

        Args:
            line (dict):  Dict containing the fields of a line.
        """
        date = "{year}-{month:02d}-{day:02d}".format(
            year=int(self.vars["yyyy"]), month=int(self.vars["m"]), day=int(self.vars["d"])
        )  # TODO: What would be a better solution for getting current date?

        for k, v in line.items():
            self.meta.setdefault(date, dict()).update({k: v})

    def _parse_velocity(self, line, cache):
        """Parse orbit velocity

        Be aware: Is not included. Needs updating to new dataset.

        Args:
            line (dict):  Dict containing the fields of a line.
            cache (dict): Temporary dictionary with the fields 'key' and 'values'.
        """
        return  # TODO: Has to be implemented correctly!!!!

        self.data.setdefault("sat_vel", list()).append(
            np.array([float(line["vel_x"]), float(line["vel_y"]), float(line["vel_z"])]) * unit.decimeter2meter
        )
        self.data.setdefault("sat_clock_rate", list()).append(float(line["clk_rate"]) * unit.microsecond2second)

        # Check if sigmas are defined in SP3 file (TODO: better solution?)
        if line["sig_vel_x"] != "":
            self.data.setdefault("sat_vel_sigma", list()).append(
                np.array([float(line["sig_vel_x"]), float(line["sig_vel_y"]), float(line["sig_vel_z"])])
                * 10 ** -4
                * self.meta[date]["base_posvel"]
                * unit.millimeter2meter
            )
            self.data.setdefault("sat_clock_rate_sigma", list()).append(
                float(line["sig_clk_rate"]) * 10 ** -4 * self.meta[date]["base_clkrate"] * unit.picosecond2second
            )

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return []

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Write data based on GNSS SP3 orbit file

        TODO:
            Add 'vel' and 'sat_clock_rate' entries to Dataset.

        Args:
            dset (Dataset): Dataset with following fields:

        ====================  ===============  =======  ========================================================
         Field                 Type             Unit     Description
        ====================  ===============  =======  ========================================================
         sat_clock_bias        numpy.ndarray     m       Satellite clock offset from GPS time
         sat_pos               PositionTable     m       Satellite position
         satellite             numpy.ndarray             Satellite PRN number
         system                numpy.ndarray             GNSS identifier
         time                  TimeTable                 Observation epoch in GPS time
        ====================  ===============  =======  ========================================================

            and following Dataset `meta` data:

        ====================  ===========================================================================
         Entry                 Description
        ====================  ===========================================================================
        agency                 Agency responsible for generating the SP3 file
        base_clkrate           Base number used for computing standard deviation of clock bias and clock rate
        base_posvel            Base number used for computing standard deviation of position and velocity
        coord_sys              Coordinate system
        data_used              Data used
        day                    Day of Gregorian date of first orbit epoch
        epoch_interval         Epoch interval between observation entries
        file_type              File type (G - GPS only, M - mixed files, R - GLONASS only, L - LEO, E - GALILEO)
        gpssec                 GPS seconds (seconds of GPS week) at the first orbit epoch
        gpsweek                GPS week at the first orbit epoch
        hour                   Hour of Gregorian date of first orbit epoch
        minute                 Minute of Gregorian date of first orbit epoch
        mjd_frac               Fractional part of Modified Julian Day at the first orbit epoch
        mjd_int                Integer part of Modified Julian Day at the first orbit epoch
        month                  Month of Gregorian date of first orbit epoch
        num_epoch              Number of epochs in the ephermeris file
        orb_type               Orbit type  (F - fitted, E - extrapolated or predicted, B - broadcast, HLM - Helmert ...)
        pv_flag                Position (P) or velocity (V) flag
        second                 Seconds of Gregorian date of first orbit epoch
        time_sys               Time system (GPS, GLO, GAL, TAI, UTC)
        version                Format version (e.g. a (GPS only), b (GPS, GLONASS), c, d)
        year                   Year of Gregorian date of first orbit epoch
        ====================  ===========================================================================

        """
        dset.num_obs = len(self.data["time"])
        dset.meta.update(self.meta)
        date = "{year}-{month:02d}-{day:02d}".format(
            year=int(self.vars["yyyy"]), month=int(self.vars["m"]), day=int(self.vars["d"])
        )  # TODO: What would be a better solution for getting current date?
        if dset.meta[date]["time_sys"] == "GPS":
            dset.add_time("time", val=self.data["time"], scale="gps")
        elif dset.meta[date]["time_sys"] == "UTC":
            dset.add_time("time", val=self.data["time"], scale="utc")
        else:
            log.fatal("Time system {} is not handled so far in Where.", dset.meta[date]["time_sys"])

        dset.add_text("satellite", val=self.data["satellite"])
        dset.add_text("system", val=self.data["system"])
        dset.add_position("sat_pos", time="time", itrs=np.array(self.data["sat_pos"]), unit="meter")
        dset.add_float("sat_clock_bias", val=np.array(self.data["sat_clock_bias"]), unit="meter")
