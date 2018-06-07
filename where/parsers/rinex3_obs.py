"""A parser for reading RINEX format 3.03 data

Example:
--------

    from where import parsers
    parser = parsers.parse('rinex3_obs', rundate=rundate, station=station)

Description:
------------

Reads data from files in the RINEX file format version 3.03 (see :cite:`rinex3`).


"""

# Standard library imports
from datetime import timedelta
import dateutil.parser
import itertools

# External library imports
import numpy as np

# Where imports
from where.parsers import parser
from where.lib import config
from where.lib import log
from where.lib import plugins
from where.lib.unit import unit


SYSTEM_TIME_OFFSET_TO_GPS_TIME = dict(BDT=14, GAL=0, IRN=0, QZS=0)


@plugins.register
class Rinex3Parser(parser.ParserDict):
    """A parser for reading RINEX observation file

    The parser reads GNSS observations in RINEX format 3.03 (see :cite:`rinex3`). The GNSS observations
    are sampled after sampling rate definition in configuration file.

    Attributes:
        data:           Dict containing the (observation) data read from file.
        dset:           Dataset object, which includes observation type (e.g. L1, C1, P1, ...), position ('site_pos')
                        and time ('time') fields.
        file_key:       Key to the RINEX observation file defined in files.conf file.
        meta:           Dict containing the metainformation read from file.
        obstypes_all:   List with all observations types for all GNSS
        sampling_rate:  Observation sampling rate in [seconds] given in configuration file.
        vars:           Variables needed to identify RINEX observation file based on definition in files.conf file.
    """

    def __init__(self, rundate=None, station=None, file_path=None, sampling_rate=None):
        """Initialize Rinex3-parser

        The file to be parsed can be specified in two different ways:

        1) For a regular Where-model run the `rundate` and `station` should be specified. The actual file will then be
           located in the archive using the files.conf configuration file.
        2) For ad-hoc parsing of a rinex-file it is also possible to give the explicit file path using the
           `file_path`-parameter.

        The sampling rate will normally be read from the configuration files, but can also be set explicitly using the
        `sampling_rate`-parameter.

        Args:
            rundate (date):           The model run date.
            station (str):            Station name (4 character station code).
            file_path (str):          Optional path to rinex-file to parse.
            sampling_rate (float):    Optional sampling rate (Default value taken from config file).
            time_scale(str):          Time scale, which is used to define the time scale of Dataset (Default: gps).
        """
        super().__init__(rundate=rundate, file_path=file_path)  # same as Parser(rundate, file_path)
        self.file_key = "gnss_rinex_obs"
        self.obstypes_all = list()
        self.time_scale = "gps"

        if station:
            self.vars["station"] = station
            self.vars["STATION"] = station.upper()

        # Sampling rate
        self.sampling_rate = config.tech.get("sampling_rate", sampling_rate).float

    #
    # PARSERS
    #
    def setup_parsers(self):
        """Parsers defined for reading RINEX observation file line by line.

           First the RINEX header information are read and afterwards the RINEX observation.
        """
        # Parser for RINEX header
        header_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line[60:73] == "END OF HEADER",
            label=lambda line, _ln: line[60:].strip(),
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      3.02           OBSERVATION DATA    M (MIXED)           RINEX VERSION / TYPE
                "RINEX VERSION / TYPE": {
                    "parser": self.parse_string,
                    "fields": {"version": (0, 20), "file_type": (20, 21), "sat_sys": (40, 41)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # MAKERINEX 2.0.20023 BKG/GOWETTZELL      2016-03-02 00:20    PGM / RUN BY / DATE
                "PGM / RUN BY / DATE": {
                    "parser": self.parse_string,
                    "fields": {"program": (0, 20), "run_by": (20, 40), "file_created": (40, 60)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # G = GPS R = GLONASS E = GALILEO S = GEO M = MIXED           COMMENT
                "COMMENT": {"parser": self.parse_comment, "fields": {"comment": (0, 60)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # stas                                                        MARKER NAME
                "MARKER NAME": {"parser": self.parse_string, "fields": {"marker_name": (0, 60)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # 66008M005                                                   MARKER NUMBER
                "MARKER NUMBER": {"parser": self.parse_string, "fields": {"marker_number": (0, 20)}},
                "MARKER TYPE": {"parser": self.parse_string, "fields": {"marker_type": (0, 20)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # SATREF              Norwegian Mapping Authority             OBSERVER / AGENCY
                "OBSERVER / AGENCY": {
                    "parser": self.parse_string, "fields": {"observer": (0, 20), "agency": (20, 60)}
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # 3008040             SEPT POLARX4        2.9.0               REC # / TYPE / VERS
                "REC # / TYPE / VERS": {
                    "parser": self.parse_string,
                    "fields": {"receiver_number": (0, 20), "receiver_type": (20, 40), "receiver_version": (40, 60)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # CR620012101         ASH701945C_M    SCIS                    ANT # / TYPE
                "ANT # / TYPE": {
                    "parser": self.parse_string, "fields": {"antenna_number": (0, 20), "antenna_type": (20, 40)}
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #   3275756.7623   321111.1395  5445046.6477                  APPROX POSITION XYZ
                "APPROX POSITION XYZ": {
                    "parser": self.parse_approx_position,
                    "fields": {"pos_x": (0, 14), "pos_y": (14, 28), "pos_z": (28, 42)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #         0.0000        0.0000        0.0000                  ANTENNA: DELTA H/E/N
                "ANTENNA: DELTA H/E/N": {
                    "parser": self.parse_float,
                    "fields": {"antenna_height": (0, 14), "antenna_east": (14, 28), "antenna_north": (28, 42)},
                },
                "ANTENNA: DELTA X/Y/Z": {
                    "parser": self.parse_float,
                    "fields": {"ant_vehicle_x": (0, 14), "ant_vehicle_y": (14, 28), "ant_vehicle_z": (28, 42)},
                },
                # TODO: 'ANTENNA:PHASECENTER'
                # TODO: 'ANTENNA:B.SIGHT XYZ'
                # TODO: 'ANTENNA:ZERODIR AZI'
                # TODO: 'ANTENNA:ZERODIR XYZ'
                # TODO: 'CENTER OF MASS: XYZ'
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # G   26 C1C C1P L1C L1P D1C D1P S1C S1P C2P C2W C2S C2L C2X  SYS / # / OBS TYPES
                #        L2P L2W L2S L2L L2X D2P D2W D2S D2L D2X S2P S2W S2S  SYS / # / OBS TYPES
                # R   16 C1C C1P L1C L1P D1C D1P S1C S1P C2C C2P L2C L2P D2C  SYS / # / OBS TYPES
                #        D2P S2C S2P                                          SYS / # / OBS TYPES
                "SYS / # / OBS TYPES": {
                    "parser": self.parse_sys_obs_types,
                    "fields": {
                        "satellite_sys": (0, 1),
                        "num_obstypes": (3, 6),
                        "type_01": (7, 10),
                        "type_02": (11, 14),
                        "type_03": (15, 18),
                        "type_04": (19, 22),
                        "type_05": (23, 26),
                        "type_06": (27, 30),
                        "type_07": (31, 34),
                        "type_08": (35, 38),
                        "type_09": (39, 42),
                        "type_10": (43, 46),
                        "type_11": (47, 50),
                        "type_12": (51, 54),
                        "type_13": (55, 58),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # DBHZ                                                        SIGNAL STRENGTH UNIT
                "SIGNAL STRENGTH UNIT": {"parser": self.parse_string, "fields": {"signal_strength_unit": (0, 20)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     1.000                                                  INTERVAL
                "INTERVAL": {"parser": self.parse_float, "fields": {"interval": (0, 10)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #  2016    03    01    00    00   00.0000000     GPS         TIME OF FIRST OBS
                "TIME OF FIRST OBS": {
                    "parser": self.parse_time_of_first_obs,
                    "fields": {
                        "year": (0, 6),
                        "month": (6, 12),
                        "day": (12, 18),
                        "hour": (18, 24),
                        "minute": (24, 30),
                        "second": (30, 43),
                        "time_sys": (48, 51),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #   2016    03    01    23    59   59.0000000     GPS         TIME OF LAST OBS
                "TIME OF LAST OBS": {
                    "parser": self.parse_time_of_last_obs,
                    "fields": {
                        "year": (0, 6),
                        "month": (6, 12),
                        "day": (12, 18),
                        "hour": (18, 24),
                        "minute": (24, 30),
                        "second": (30, 43),
                        "time_sys": (48, 51),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      0                                                      RCV CLOCK OFFS APPL
                "RCV CLOCK OFFS APPL": {"parser": self.parse_string, "fields": {"rcv_clk_offset_flag": (0, 6)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # G APPL_DCB          xyz.uvw.abc//pub/dcb_gps.dat            SYS / DCBS APPLIED
                "SYS / DCBS APPLIED": {
                    "parser": self.parse_sys_dcbs_applied,
                    "fields": {"sat_sys": (0, 1), "program": (2, 19), "source": (20, 60)},
                },
                "SYS / PCVS APPLIED": {
                    "parser": self.parse_sys_pcvs_applied,
                    "fields": {"sat_sys": (0, 1), "program": (2, 19), "source": (20, 60)},
                },
                "SYS / SCALE FACTOR": {
                    "parser": self.parse_scale_factor,  # TODO: not implemented
                    "fields": {"sat_sys": (0, 1), "factor": (2, 6), "num_obstypes": (8, 10)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # G L1C  0.00000  12 G01 G02 G03 G04 G05 G06 G07 G08 G09 G10  SYS / PHASE SHIFT
                #                    G11 G12                                  SYS / PHASE SHIFT
                # G L1W  0.00000                                              SYS / PHASE SHIFT
                "SYS / PHASE SHIFT": {
                    "parser": self.parse_phase_shift,
                    "fields": {
                        "sat_sys": (0, 1),
                        "obs_type": (2, 5),
                        "correction": (6, 14),
                        "num_satellite": (16, 18),
                        "satellites": (19, 59),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #  22 R01  1 R02 -4 R03  5 R04  6 R05  1 R06 -4 R07  5 R08  6 GLONASS SLOT / FRQ #
                #     R09 -6 R10 -7 R11  0 R13 -2 R14 -7 R15  0 R17  4 R18 -3 GLONASS SLOT / FRQ #
                #     R19  3 R20  2 R21  4 R22 -3 R23  3 R24  2               GLONASS SLOT / FRQ #
                "GLONASS SLOT / FRQ #": {
                    "parser": self.parse_glonass_slot,
                    "fields": {
                        "num_satellite": (0, 3),
                        "slot_01": (4, 11),
                        "slot_02": (11, 18),
                        "slot_03": (18, 25),
                        "slot_04": (25, 32),
                        "slot_05": (32, 39),
                        "slot_06": (39, 46),
                        "slot_07": (46, 53),
                        "slot_08": (53, 60),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #  C1C  -10.000 C1P  -10.123 C2C  -10.432 C2P  -10.634        GLONASS COD/PHS/BIS
                "GLONASS COD/PHS/BIS": {
                    "parser": self.parse_glonass_code_phase_bias,
                    "fields": {"type_01": (1, 13), "type_02": (14, 26), "type_03": (27, 39), "type_04": (40, 52)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     16    17  1851     3                                    LEAP SECONDS
                "LEAP SECONDS": {
                    "parser": self.parse_leap_seconds,
                    "fields": {
                        "leap_seconds": (0, 6),
                        "future_past_leap_seconds": (6, 12),
                        "week": (12, 18),
                        "week_day": (18, 24),
                        "time_sys": (24, 27),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     71                                                      # OF SATELLITES
                "# OF SATELLITES": {"parser": self.parse_integer, "fields": {"num_satellites": (0, 6)}},
                # TODO: 'PRN / # OF OBS'
            },
        )

        # Parser for RINEX observation blocks

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9
        # > 2006 03 24 13 10 36.0000000  0  5      -0.123456789012
        # G06  23629347.915            .300 8         -.353 4  23629347.158          24.158
        # G09  20891534.648           -.120 9         -.358 6  20891545.292          38.123
        # E11          .324 8          .178 7
        # S20  38137559.506      335849.135 9
        obs_parser = parser.define_parser(
            end_marker=lambda _l, _ln, next_line: next_line.startswith(">"),
            label=lambda line, _ln: (
                line[0:1].isalpha()  # Obs line start with sat. system identifier
                and not line.startswith(">")  # Observation epoch line starts with '>'
                and not line[60:61].isalpha()
            ),  # Comment line
            parser_def={
                False: {
                    "parser": self.parse_observation_epoch,
                    "fields": {
                        "year": (2, 6),
                        "month": (7, 9),
                        "day": (10, 12),
                        "hour": (13, 15),
                        "minute": (16, 18),
                        "second": (18, 29),
                        "epoch_flag": (31, 32),
                        "num_sat": (32, 35),
                        "rcv_clk_offset": (41, 56),
                        "comment": (60, 80),
                    },
                },
                True: {
                    "parser": self.parse_observation,
                    "strip": "\n",  # Remove only newline '\n' leading and trailing characters from line.
                    "fields": {
                        "sat": (0, 3), "obs": (3, None)  # 'None' indicates, that line is sliced until the end of line.
                    },
                },
            },
        )

        return itertools.chain([header_parser], itertools.repeat(obs_parser))

    #
    # HEADER PARSERS
    #
    def parse_approx_position(self, line, _):
        """Parse station coordinates defined in RINEX header to instance variable `data`.
        """
        pos = np.array((float(line["pos_x"]), float(line["pos_y"]), float(line["pos_z"])))
        self.data["pos"] = pos

    def parse_comment(self, line, _):
        """Parse comment lines in RINEX header to instance variable `meta['comment']`.
        """
        self.meta.setdefault("comment", list()).append(line["comment"])

    def parse_float(self, line, _):
        """Parse float entries of RINEX header to instance variable `meta`.
        """
        self.parse_default_meta({k: float(v) for k, v in line.items()}, _)

    def parse_glonass_code_phase_bias(self, line, _):
        """Parse GLONASS phase correction in RINEX header to instance variable `meta['glonass_bias']`.

            self.meta['glonass_bias'] = { <obstype>: <bias in meters>}
        """
        self.meta.setdefault("glonass_bias", {})
        for field in sorted([f for f in line if f.startswith("type_")]):
            if line[field]:
                type_, bias = line[field].split()[0:2]
                self.meta["glonass_bias"].update({type_: float(bias)})

    def parse_glonass_slot(self, line, _):
        """Parse GLONASS slot and frequency numbers given in RINEX header to instance variable `meta['glonass_slot']`.

            self.meta['glonass_slot'] = { <slot>: <frequency number>}
        """
        self.meta.setdefault("glonass_slot", {})
        for field in sorted([f for f in line if f.startswith("slot_")]):
            if line[field]:
                slot, freq = line[field].split()[0:2]
                self.meta["glonass_slot"].update({slot: int(freq)})

    def parse_integer(self, line, _):
        """Parse integer entries of RINEX header to instance variable `meta`.
        """
        self.parse_default_meta({k: int(v) for k, v in line.items()}, _)

    def parse_leap_seconds(self, line, _):
        """Parse entries of RINEX header `LEAP SECONDS` to instance variable `meta`.

            self.meta['leap_seconds'] = { 'leap_seconds': <value>,
                                          'future_past_leap_seconds': <value>,
                                          'week': <value>,
                                          'week_day': <value>,
                                          'time_sys': <system> }
        """
        for field in line:
            self.meta.setdefault("leap_seconds", {}).update({field: line[field]})

    def parse_phase_shift(self, line, cache):
        """Parse entries of RINEX header `SYS / PHASE SHIFT` to instance variable `meta`.

            self.meta['phase_shift'] = { <sat_sys>: { <obs_type>: { corr: <correction>,
                                                                    sat: <[satellite list]>}}}

        Example of `phase_shift` meta entry:

            self.meta['phase_shift'] =  {'G': {'L1C': {'corr': '0.00000',
                                                       'sat': ['G01', 'G02', 'G03', ...]},
                                               'L1W': {'corr': '0.00000',
                                                       'sat': []}},
                                         'R': {'L1C': {'corr': '0.00000',
                                                       'sat': ['R01', 'R02', 'R07', 'R08']}}}

        TODO: Maybe better to add information to meta['obstypes']?
        """
        self.meta.setdefault("phase_shift", {})

        if line["sat_sys"]:
            cache["sat_sys"] = line["sat_sys"]
            cache["obs_type"] = line["obs_type"]
            cache["corr"] = line["correction"]
            cache["sat"] = []

        if cache["sat_sys"] not in self.meta["phase_shift"]:
            self.meta["phase_shift"].update({cache["sat_sys"]: {}})

        cache["sat"].extend(line["satellites"].split())

        if cache["obs_type"]:
            self.meta["phase_shift"][cache["sat_sys"]].update(
                {cache["obs_type"]: {"corr": cache["corr"], "sat": cache["sat"]}}
            )

    def parse_scale_factor(self, line, _):
        """Parse entries of RINEX header `SYS / SCALE FACTOR` to instance variable `meta`.
        """
        log.fatal("Reading and applying of RINEX header entry 'SYS / SCALE FACTOR' is not implemented.")

    def parse_string(self, line, _):
        """Parse string entries of RINEX header to instance variable 'meta'.
        """
        self.parse_default_meta({k: v for k, v in line.items()}, _)

    def parse_sys_dcbs_applied(self, line, _):
        """Parse entries of RINEX header `SYS / DCBS APPLIED` to instance variable `meta`.

            self.meta['dcbs_applied'] = { <sat_sys>: { prg: <used program>,
                                                       url: <source url>}}
        """
        self.meta.setdefault("dcbs_applied", {}).update(
            {line["sat_sys"]: {"prg": line["program"], "url": line["source"]}}
        )

    def parse_sys_obs_types(self, line, cache):
        """Parse observation types given in RINEX header to instance variable `meta['obstypes']` and data.

        The data dictionaries `obs`, `cycle_slip` and `signal_strength` are initialized based on the given observation
        type in the RINEX header.

            self.meta['obstypes'] = { <sat_sys>: [<ordered list with given observation types>]}
        """
        self.data.setdefault("obs", {})
        self.data.setdefault("cycle_slip", {})
        self.data.setdefault("signal_strength", {})
        self.meta.setdefault("obstypes", {})

        if line["satellite_sys"]:
            cache["sys"] = line["satellite_sys"]
            cache["obstypes"] = list()

        for field in sorted([f for f in line if f.startswith("type_")]):
            if line[field]:
                cache["obstypes"].append(line[field])
                if line[field] not in self.obstypes_all:
                    self.obstypes_all.append(line[field])
                self.meta["obstypes"].update({cache["sys"]: cache["obstypes"]})
                self.data["obs"][line[field]] = list()
                self.data["cycle_slip"][line[field]] = list()
                self.data["signal_strength"][line[field]] = list()

    def parse_sys_pcvs_applied(self, line, _):
        """Parse entries of RINEX header `SYS / PCVS APPLIED` to instance variable `meta`.

            self.meta['pcvs_applied'] = { <sat_sys>: { prg: <used program>,
                                                       url: <source url>}}
        """
        self.meta.setdefault("pcvs_applied", {}).update(
            {line["sat_sys"]: {"prg": line["program"], "url": line["source"]}}
        )

    def parse_time_of_first_obs(self, line, _):
        """Parse time of first observation given in RINEX header to instance variable `meta`.
        """
        if line["time_sys"] != "GPS":
            log.fatal("Time system {} is not handled so far in Where.", line["time_sys"])

        if line["time_sys"] not in self.meta:
            self.meta["time_sys"] = line["time_sys"]
        else:
            if line["time_sys"] != self.meta["time_sys"]:
                log.fatal(
                    "Time system definition in 'TIME OF FIRST OBS' ({}) and 'TIME OF LAST OBS' ({}) are not"
                    "identical.",
                    line["time_sys"],
                    self.meta["time_sys"],
                )

        if line["year"]:
            self.meta["time_first_obs"] = (
                "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
                "".format(
                    year=int(line["year"]),
                    month=int(line["month"]),
                    day=int(line["day"]),
                    hour=int(line["hour"]),
                    minute=int(line["minute"]),
                    second=float(line["second"]),
                )
            )

    def parse_time_of_last_obs(self, line, _):
        """Parse time of last observation given in RINEX header to instance variable `meta`.
        """
        if line["time_sys"] != "GPS":
            log.fatal("Time system {} is not handled so far in Where.", line["time_sys"])

        if line["time_sys"]:
            if line["time_sys"] not in self.meta:
                self.meta["time_sys"] = line["time_sys"]
            else:
                if line["time_sys"] != self.meta["time_sys"]:
                    log.fatal(
                        "Time system definition in 'TIME OF FIRST OBS' ({}) and 'TIME OF LAST OBS' ({}) are"
                        "not identical.",
                        self.meta["time_sys"],
                        line["time_sys"],
                    )

        if line["year"]:
            self.meta["time_last_obs"] = (
                "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
                "".format(
                    year=int(line["year"]),
                    month=int(line["month"]),
                    day=int(line["day"]),
                    hour=int(line["hour"]),
                    minute=int(line["minute"]),
                    second=float(line["second"]),
                )
            )

    #
    # OBSERVATION PARSERS
    #
    def parse_observation_epoch(self, line, cache):
        """Parse observation epoch information of RINEX observation record

        In addition the RINEX observation are decimated based on the given sampling rate.
        """
        # Reject empty line -> TODO: Better solution?
        if not line["year"].isnumeric():
            return

        # Reject comment line
        if line["comment"][0:1].isalpha():
            return

        cache["obs_time"] = (
            "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
            "".format(
                year=int(line["year"]),
                month=int(line["month"]),
                day=int(line["day"]),
                hour=int(line["hour"]),
                minute=int(line["minute"]),
                second=float(line["second"]),
            )
        )
        cache["obs_sec"] = (
            int(line["hour"]) * unit.hour2second + int(line["minute"]) * unit.minute2second + float(line["second"])
        )
        cache["epoch_flag"] = int(line["epoch_flag"])
        cache["rcv_clk_offset"] = _float(line["rcv_clk_offset"])

        if line["epoch_flag"].strip() != "0":
            log.fatal(
                "Epoch {} is not ok, which is indicated by epoch flag {}. How it should be handled in Where?",
                cache["obs_time"],
                line["epoch_flag"],
            )  # TODO: Handle flagged epochs

        # Decimate RINEX observation defined by sampling rate [seconds]
        if cache["obs_sec"] % self.sampling_rate != 0:
            cache["obs_sec"] = None  # Ignore epoch

    def parse_observation(self, line, cache):
        """Parse observation record of RINEX file

        In addition the RINEX observation are rejected from satellites not given in configuration file.
        """
        # Ignore epochs based on sampling rate
        # TODO: This should be done in 'edit' step!!!
        sec = cache["obs_sec"]
        if sec is None:
            return

        # Fit length of observation line against definition (16 characters * number of observation types)
        #
        # NOTE: This is necessary, because missing observations are written as 0.0 or BLANK in RINEX format. Therefore
        #       the last observation can be BLANK. In this case the length of observation line does not fit the
        #       definition (16 characters * number of observation types), which is assumed by reading the observation
        #       line in the following.
        sys = line["sat"][0]
        num_obstypes = len(self.meta["obstypes"][sys])
        field_length = 16
        line_length = field_length * num_obstypes
        line["obs"] = line["obs"].ljust(line_length)

        # Parse observation line in fields
        for idx, obs_type in zip(range(0, line_length, field_length), self.meta["obstypes"][sys]):
            value = line["obs"][idx:idx + field_length]
            self.data["obs"][obs_type].append(_float(value[0:14]))
            self.data["cycle_slip"][obs_type].append(_int(value[14:15]))
            self.data["signal_strength"][obs_type].append(_int(value[15:16]))

        # Fill unused observation types with zero values
        #
        # NOTE: Each observation type is saved in a Dataset field. The observation type fields have the same length
        #       to be consistent with the time, system or satellite Dataset field. The problem is that some observation
        #       types are not observed for a certain satellite system, but these observation are included with zero
        #       values in the observation type field, which is done in the following.
        unused_obstypes = set(self.obstypes_all) - set(self.meta["obstypes"][sys])
        for obs_type in unused_obstypes:
            self.data["obs"][obs_type].append(0.0)
            self.data["cycle_slip"][obs_type].append(0)
            self.data["signal_strength"][obs_type].append(0)

        self.data.setdefault("time", list()).append(cache["obs_time"])
        self.data.setdefault("epoch_flag", list()).append(cache["epoch_flag"])
        self.data.setdefault("rcv_clk_offset", list()).append(cache["rcv_clk_offset"])

        obs = {
            "station": self.meta["marker_name"].lower(),  # vars['station'],
            "site_id": self.meta["marker_name"].lower(),
            "system": sys,
            "satellite": line["sat"],
            "satnum": line["sat"][1:3],
        }

        for field, value in obs.items():
            self.data.setdefault("text", dict()).setdefault(field, list()).append(value)

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return [self.remove_empty_systems, self.remove_empty_obstype_fields, self.time_system_correction]

    def remove_empty_systems(self):
        """Remove GNSSs without observations from `self.meta['obstypes']`.

        The GNSSs are defined in RINEX header (SYS / # / OBS TYPES ). It can happen, that no observations are available
        for a given GNSS. GNSSs without observations are deleted from dictionary `self.meta['obstypes']` in this
        routine.
        """
        for sys in list(self.meta["obstypes"].keys()):
            if sys not in self.data["text"]["system"]:
                log.debug("No observation given for GNSS '{}'. GNSS '{}' is removed from Dataset.", sys, sys)
                del self.meta["obstypes"][sys]

    def remove_empty_obstype_fields(self):
        """Remove empty observation type data fields.

        The observation types given in RINEX header (SYS / # / OBS TYPES ) define the key entries initialized in the
        data dictionaries `obs`, `cycle_slip` and `signal_strength`. It can happen, that no observations are available
        for observation types defined in the RINEX header. Empty observation type entries are deleted in this routine.

        In addition it can be that the observation types are not given for a certain GNSS. In this case the observation
        types in the meta['obstypes'] variable are removed for this GNSS. But it should be kept in mind, that the
        observations for this observation type are still given in the Dataset, which are set to zero.
        """
        remove_obstype = []  # List with observation types, which should be removed from Dataset.
        remove_obstype_sys = {}  # Dictionary with observation types given for each GNSS, which should be removed from
        # meta['obstypes'].
        for obstype, obs in self.data["obs"].items():
            if not obs or np.all(np.array(obs) == 0.0):
                remove_obstype.append(obstype)
            systems = set(self.data["text"]["system"])
            for sys in systems:

                # Filter observations depending on GNSS
                idx = np.array(self.data["text"]["system"]) == sys

                if np.all(np.array(obs)[idx] == 0.0):
                    remove_obstype_sys.setdefault(sys, list()).append(obstype)

        log.debug(
            "Following observation types are removed, because no observation given: {}",
            " ".join(sorted(remove_obstype)),
        )

        # Remove empty observation type data fields
        for obstype in remove_obstype:
            for sys in list(self.meta["obstypes"]):
                if obstype in self.meta["obstypes"][sys]:
                    self.meta["obstypes"][sys].remove(obstype)

            del self.data["obs"][obstype]
            del self.data["cycle_slip"][obstype]
            del self.data["signal_strength"][obstype]

        # Remove empty observation types for a given GNSS from meta['obstypes'] and other meta variables
        for sys, obstypes in remove_obstype_sys.items():
            for obstype in obstypes:
                if obstype in self.meta["obstypes"][sys]:
                    self.meta["obstypes"][sys].remove(obstype)

    def time_system_correction(self):
        """Apply correction to given time system for getting GPS or UTC time scale

        Following relationship are given between GNSS time scale (either BeiDou, Galileo, IRNSS or QZSS)
        :math:`t_{GNSS}` and GPS time scale :math:`t_{GPS}` (see Section 2.1.4 in :cite:`teunissen2017`):
        .. math::
              t_{GPS}  = t_{GNSS} + \Delta t

        The time offset :math:`\Delta t` is 0 s for Galileo, IRNSS and QZSS and for BeiDou 14 s. All these time scales
        are related to the International Atomic Time (TAI) by a certain time offset. An exception is the GLONASS time
        scale, which is related to UTC:
        .. math::
              t_{UTC}  = t_{GLONASS} - 3h

        Note, that in the RINEX format (see section 8.1 in :cite:`rinex3`) GLONASS time has the same hours as UTC and
        not UTC + 3h as the original GLONASS system time, which is given in the Moscow time zone instead of Greenwich.

        In this routine the given observation time (epoch) will be transformed to GPS time scale for BeiDou, Galileo,
        QZSS and IRNSS and to UTC time scale for GLONASS.
        """
        system = self.meta["time_sys"]
        valid_time_systems = ["BDT", "GAL", "GPS", "GLO", "IRN", "QZS"]

        if system not in valid_time_systems:
            log.fatal(
                "Time system '{}' in file {} is not handled in Where. Following time systems can be used: {}.",
                system,
                self.file_path,
                ", ".join(valid_time_systems),
            )

        # Convert observation time entries of BeiDou to GPS time scale by adding system time offset
        if system == "BDT":
            self.data["time"] = [
                dateutil.parser.parse(t) + timedelta(seconds=SYSTEM_TIME_OFFSET_TO_GPS_TIME.get(system, 0))
                for t in self.data["time"]
            ]

        # Change time scale to UTC for GLONASS
        elif system == "GLO":
            self.time_scale = "utc"

    # def pseudorange_system_correction(self):
    #    """Apply correction to pseudorange observations
    #
    #    See section 8.2 in :cite:`rinex3`.
    #    TODO: Is it necessary to correct the pseudorange observation depending on used time system?
    #    """

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Store GNSS data in a dataset

        Args:
            dset (Dataset): The Dataset where GNSS observation are stored with following fields:

        ====================  ==================  =================================================================
         Field                 Type                Description
        ====================  ==================  =================================================================
        <observation type>     numpy.ndarray       GNSS observation type data (e.g. C1C, C2W, L1C, L2W, ...) given
                                                   in meters
        epoch_flag             numpy.ndarray       Epoch flag
        rcv_clk_offset         numpy.ndarray       Receiver clock offset in seconds given for each epoch
        satellite              numpy.ndarray       Satellite PRN number together with GNSS identifier (e.g. G07)
        satnum                 numpy.ndarray       Satellite PRN number (e.g. 07)
        site_pos               PositionTable       PositionTable object with given station coordinates (read from
                                                   RINEX header)
        station                numpy.ndarray       Station name list
        system                 numpy.ndarray       GNSS identifier
        time                   TimeTable           Observation time given as TimeTable object
        ====================  ==================  =================================================================

            and following Dataset `meta` data:

        ==================== ======= ===============================================================================
         Entry                Type    Description
        ==================== ======= ===============================================================================
        agency               str     Name of agency from observer
        antenna_east         float   East component of vector between marker and antenna reference point in meters
        antenna_height       float   Height component of vector between marker and antenna reference point in meters
        antenna_north        float   North component of vector between marker and antenna reference point in meters
        antenna_number       str     Antenna serial number
        antenna_type         str     Antenna type
        ant_vehicle_x        float   X-coordinate in body-fixed coord. system of antenna reference point on vehicle
        ant_vehicle_y        float   Y-coordinate in body-fixed coord. system of antenna reference point on vehicle
        ant_vehicle_z        float   Z-coordinate in body-fixed coord. system of antenna reference point on vehicle
        comment              list    List with RINEX header comment lines
        dcbs_applied         dict    Satellite system dependent information about applying DCBs
        file_created         str     Date and time of file creation
        file_type            str     File type (e.g. 'O' for observation data)
        glonass_bias         dict    GLONASS phase bias correction in meters given for code observation type (C1C, C1P,
                                     C2C and/or C2P)
        glonass_slot         dict    GLONASS frequency numbers given for GLONASS slot
        interval             float   Observation interval in seconds
        leap_seconds         dict    Dictionary with information related to leap seconds
        marker_name          str     Name of antenna marker
        marker_number        str     Number of antenna marker
        num_satellites       int     Number of satellites, for which observations are stored in the RINEX file
        observer             str     Name of observer
        obstypes             dict    Observation types given for each GNSS
        pcvs_applied         dict    Satellite system dependent information about applying PCVs
        phase_shift          dict    Phase shift correction given for a satellite system dependent observation type
        program              str     Name of program creating current file
        rcv_clk_offset_flag  str     Flag (1=yes, 0=no) indicating if realtime-derived receiver clock offset is
                                     applied for epoch, code, and phase
        receiver_number      str     Receiver serial number
        receiver_type        str     Receiver type
        receiver_version     str     Receiver firmware version
        run_by               str     Name of agency creating current file
        sat_sys              str     Satellite system given in observation file (G, R, E, J, C, I, S or M)
        signal_strength_unit str     Unit of the carrier to noise ratio observables
        time_first_obs       str     Time of first observation record
        time_last_obs        str     Time of last observation record
        time_sys             str     Time system used for GNSS observations (GPS, GLO, GAL, QZS, BDT or IRN)
        version              str     Format version
        ==================== ======= ===============================================================================
        """
        # Time
        dset.num_obs = len(self.data["time"])
        dset.add_time("time", val=self.data["time"], scale=self.time_scale)
        dset.add_float("epoch_flag", val=np.array(self.data["epoch_flag"]))
        dset.add_float("rcv_clk_offset", val=np.array(self.data["rcv_clk_offset"]))

        dset.meta.update(self.meta)

        for field, value in self.data["text"].items():
            dset.add_text(field, val=value)

        # Observations
        for obs_type in self.data["obs"]:
            dset.add_float(obs_type, table="obs", val=np.array(self.data["obs"][obs_type]), unit="meter")
            dset.add_float(obs_type + "_lli", table="lli", val=np.array(self.data["cycle_slip"][obs_type]))
            dset.add_float(obs_type + "_snr", table="snr", val=np.array(self.data["signal_strength"][obs_type]))

        # Positions
        dset.add_position("site_pos", time="time", itrs=np.repeat(self.data["pos"][None, :], dset.num_obs, axis=0))


def _float(value):
    """Convert string to float value

    Whitespace or empty value is set to 0.0.

    Args:
        value (str):   string value

    Returns:
        float: float value
    """
    if value.isspace() or not value:
        return 0.0
    else:
        return float(value)


def _int(value):
    """Convert string to int value

    Whitespace or empty value is set to 0.

    Args:
        value (str): string value

    Returns:
        int: integer value
    """
    if value.isspace() or not value:
        return 0
    else:
        return int(value)
