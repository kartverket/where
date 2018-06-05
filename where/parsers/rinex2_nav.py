"""A parser for reading GPS RINEX navigation file (exception GLONASS and SBAS)

Example:
--------

    from where import data
    from where import parsers

    # Parse data
    parser = parsers.parse('rinex2_nav', rundate=rundate)

    # Create a empty Dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=0, empty=True)

    # Fill Dataset with parsed data
    parser.write_to_dataset(dset)

Description:
------------

Reads GNSS data from files in the RINEX navigation file format 2.11 (see :cite:`rinex2`). Note that beside the
given day also the navigation files from the day before and the day after is read. An exception is also, that this
parser does not handle GLONASS and SBAS navigation messages. All navigation time epochs (time of clock (toc)) are
converted to GPS time scale.

The navigation message is not defined for GALILEO, BeiDou, QZSS and IRNSS in RINEX format 2.11. In this case the RINEX
3.03 definition is used (see :cite:`rinex3`).


$Revision: 15027 $
$Date: 2018-05-08 15:26:26 +0200 (Tue, 08 May 2018) $
$LastChangedBy: dahmic $
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
from where.lib import files
from where.lib import gnss
from where.lib import log
from where.lib import plugins
from where.lib.time import Time, TimeDelta

# The constant shows the 'second' relation between the GNSS time systems to the GPS time scale. Galileo (E),
# IRNSS (I) and QZSS (J) uses the same second as the GPS time scale, but BeiDou (C) BDT time system is 14 seconds behind
# GPS time.
SYSTEM_TIME_OFFSET_TO_GPS_SECOND = dict(C=14, E=0, I=0, J=0)

# The constant shows the relation between the GNSS week to the GPS week. Galileo (E), IRNSS (I) and QZSS (J) week
# corresponds to GPS week, whereas BeiDou (C) week starts at GPS week 1356.
SYSTEM_TIME_OFFSET_TO_GPS_WEEK = dict(C=1356, E=0, I=0, J=0)

# The constant shows the relation between the navigation file extensions and the GNSS.
SYSTEM_FILE_EXTENSION = {"n": "G", "g": "R", "l": "E"}

SYSNAMES = dict(
    gnss_data_info={"G": "codes_l2", "J": "codes_l2", "E": "data_source"},
    gnss_interval={"G": "fit_interval", "J": "fit_interval", "C": "age_of_clock_corr"},
    gnss_iodc_groupdelay={"G": "iodc", "J": "iodc", "E": "bgd_e1_e5b", "C": "tgd_b2_b3"},
    gnss_l2p_flag={"G": "l2p_flag", "J": "l2p_flag"},
    gnss_tgd_bgd={"G": "tgd", "J": "tgd", "E": "bgd_e1_e5a", "C": "tgd_b1_b3", "I": "tgd"},
)


@plugins.register
class Rinex2NavParser(parser.Parser):
    """A parser for reading RINEX navigation file

    The parser reads GNSS broadcast ephemeris in RINEX format 2.11 (see :cite:`rinex2`).

    #TODO: Would it not be better to use one leading underscore for non-public methods and instance variables.

    Attributes:
        data (dict):             Dict containing the (navigation) data read from file.
        data_available (bool):   Indicator of whether data are available.
        day_offset (int):        Day offset used to calculate the number of days to read.
        dependencies (list):     List of files that have been read by the parser.
        file_key (str):          Key to the SP3 orbit file defined in files.conf file.
        file_path (pathlib.PosixPath):  File path to SP3 orbit file.
        meta (dict):             Dict containing the metainformation read from file, whereby the current day define the
                                 key and the values are the RINEX header entries.
        rundate (datetime.date): Date of model run.
        system(str):             GNSS identifier.
        vars (dict):             Variables needed to identify RINEX observation file based on definition in files.conf
                                 file.

    Methods:
        calculate_data()          Carry out the defined calculators
        copy_cache_to_data()      Copy contents of the cache to the data datastructure
        copy_cache_to_meta()      Copy contents of the cache to the meta datastructure
        parse()                   Parse data
        parse_default()           Add the contents of line to data
        parse_default_meta()      Add the contents of line to meta
        parse_line()              Parse line
        read_data()               Read data from datafiles
        setup_calculators()       List steps necessary for postprocessing
        setup_parsers()           Setup parser definition
        write_to_dataset()        Write data based on GNSS SP3 orbit file

        _check_nav_message()      Check correctness of navigation message
        _get_system_from_file_extension()  Get GNSS by reading RINEX navigation file extension
        _parse_file()             Read a data file and parse the content
        _parse_float()            Parse float entries of SP3 header to instance variable 'meta'
        _parse_ion_alpha()        Parse entries of RINEX header `ION ALPHA` to instance variable `meta`.
        _parse_ion_beta()         Parse entries of RINEX header `ION BETA` to instance variable `meta`.
        _parse_obs_float()        Parse float entries of RINEX navigation data block to instance variable 'data'.
        _parse_observation_epoch()  Parse observation epoch information of RINEX navigation data record
        _parse_string()           Parse string entries of SP3 header to instance variable 'meta'
        _parse_string_list()      Parse string entries of RINEX header to instance variable 'meta' in a list
        _parse_time_system_corr() Not defined.
        _rename_fields_based_on_system()  Rename general GNSS fields to GNSS specific ones
        _time_system_correction() Apply correction to given time system for getting GPS time scale
    """

    def __init__(self, rundate, station):
        """
        Args:
            rundate (datetime.datetime):   The model run date.
            station (str):                 Station name.
        """
        super().__init__(rundate)  # same as Parser(rundate)
        self.vars["station"] = station
        self.vars["STATION"] = station.upper()

    def read_data(self):
        """Read RINEX navigation files

        Overridden to read data for three consecutive days, the given date and the day before.

        TODO: Is that necessary to read three consecutive days of broadcast ephemeris???? No interpolation needed. The
              valid RINEX epoch is used for given time. The routine get_brdc_idx(dset) has to be changed.
        """
        self.system = self._get_system_from_file_extension()  # TODO: Check design of solution?
        date_to_read = self.rundate - timedelta(days=1)

        while date_to_read <= (self.rundate + timedelta(days=1)):
            self.vars.update(config.date_vars(date_to_read))
            self.vars["date"] = date_to_read.strftime("%Y-%m-%d")
            file_path = files.path(self.file_key, file_vars=self.vars)
            self.dependencies.append(file_path)

            with files.open_path(file_path, mode="rt") as fid:
                self.parse_file(fid)

            date_to_read += timedelta(days=1)

    def _get_system_from_file_extension(self):
        """Get GNSS by reading RINEX navigation file extension

        For RINEX 2.11 navigation files is only the file extension an indicator, which GNSS is used.

        """
        file_extension = str(self.file_path)[-1].lower()
        try:
            return SYSTEM_FILE_EXTENSION[file_extension]
        except KeyError:
            log.fatal(
                "RINEX 2.11 navigation file extension '{}' of file {} does not exists or is not handled by "
                "Where. Following RINEX navigation file extensions can be used: {}.",
                file_extension,
                self.file_path,
                ", ".join(SYSTEM_FILE_EXTENSION.keys()),
            )

    #
    # PARSERS
    #
    def setup_parsers(self):
        """Parsers defined for reading RINEX navigation file line by line

           First the RINEX header information are read and afterwards the RINEX navigation data block.
        """
        # Parser for RINEX header
        header_parser = parser.define_parser(
            end_marker=lambda line, _ln, _n: line[60:73] == "END OF HEADER",
            label=lambda line, _ln: line[60:].strip(),
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      2              NAVIGATION DATA                         RINEX VERSION / TYPE
                "RINEX VERSION / TYPE": {
                    "parser": self._parse_string, "fields": {"version": (0, 20), "file_type": (20, 21)}
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # CCRINEXN V1.6.0 UX  CDDIS               01-MAR-16 19:39     PGM / RUN BY / DATE
                "PGM / RUN BY / DATE": {
                    "parser": self._parse_string,
                    "fields": {"program": (0, 20), "run_by": (20, 40), "file_created": (40, 60)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # IGS BROADCAST EPHEMERIS FILE                                COMMENT
                "COMMENT": {"parser": self._parse_string_list, "fields": {"comment": (0, 60)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     0.1583D-07  0.0000D+00 -0.1192D-06  0.0000D+00          ION ALPHA
                "ION ALPHA": {
                    "parser": self._parse_ion_alpha,
                    "fields": {"ion_a0": (0, 14), "ion_a1": (14, 26), "ion_a2": (26, 38), "ion_a3": (38, 50)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     0.1126D+06 -0.1638D+05 -0.2621D+06  0.6554D+05          ION BETA
                "ION BETA": {
                    "parser": self._parse_ion_beta,
                    "fields": {"ion_b0": (0, 14), "ion_b1": (14, 26), "ion_b2": (26, 38), "ion_b3": (38, 50)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #    -0.931322574615D-09-0.444089209850D-14   405504     1886 DELTA-UTC: A0,A1,T,W
                #
                "DELTA-UTC: A0,A1,T,W": {
                    "parser": self._parse_time_system_corr,
                    "fields": {"a0": (0, 22), "a1": (22, 41), "t": (41, 50), "w": (50, 59)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #   1998     2    16    0.379979610443D-06                    CORR TO SYSTEM TIME
                # TODO: Only used by GLONASS
                # TODO 'CORR TO SYSTEM TIME':
                #    {'parser': self.parse_corr_to_system_time,
                #     'fields': {'year':                (0, 6),
                #                'month'               (6, 12),
                #                'day':               (12, 18),
                #                'corr':              (21, 40),
                #               }},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     17                                                      LEAP SECONDS
                "LEAP SECONDS": {"parser": self._parse_leap_seconds, "fields": {"leap_seconds": (0, 6)}},
            },
        )

        # Parser for RINEX navigation data blocks

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        #  1 16  3  1  0  0  0.0 0.130985863507D-04 0.102318153950D-11 0.000000000000D+00
        #     0.300000000000D+02-0.152187500000D+02 0.451304504878D-08-0.111734822451D+01
        #    -0.806525349617D-06 0.509567500558D-02 0.698119401932D-05 0.515365027618D+04
        #     0.172800000000D+06 0.119209289551D-06-0.654720735518D+00-0.540167093277D-07
        #     0.963950790865D+00 0.247875000000D+03 0.460218159465D+00-0.807676503456D-08
        #     0.664313395959D-10 0.100000000000D+01 0.188600000000D+04 0.000000000000D+00
        #     0.200000000000D+01 0.000000000000D+00 0.512227416039D-08 0.300000000000D+02
        #     0.172800000000D+06 0.000000000000D+00 0.000000000000D+00 0.000000000000D+00
        data_parser = parser.define_parser(
            end_marker=lambda _l, line_num, _n: line_num == 8,
            label=lambda _l, line_num: line_num,
            parser_def={
                1: {
                    "parser": self._parse_observation_epoch,
                    "fields": {
                        "sat": (0, 2),
                        "year": (2, 5),
                        "month": (5, 8),
                        "day": (8, 11),
                        "hour": (11, 14),
                        "minute": (14, 17),
                        "second": (17, 22),
                        "sat_clock_bias": (22, 41),
                        "sat_clock_drift": (41, 60),
                        "sat_clock_drift_rate": (60, 79),
                    },
                },
                2: {
                    "parser": self._parse_obs_float,
                    "fields": {"iode": (0, 22), "crs": (22, 41), "delta_n": (41, 60), "m0": (60, 79)},
                },
                3: {
                    "parser": self._parse_obs_float,
                    "fields": {"cuc": (0, 22), "e": (22, 41), "cus": (41, 60), "sqrt_a": (60, 79)},
                },
                4: {
                    "parser": self._parse_obs_float,
                    "fields": {"toe": (0, 22), "cic": (22, 41), "Omega": (41, 60), "cis": (60, 79)},
                },
                5: {
                    "parser": self._parse_obs_float,
                    "fields": {"i0": (0, 22), "crc": (22, 41), "omega": (41, 60), "Omega_dot": (60, 79)},
                },
                6: {
                    "parser": self._parse_obs_float,
                    "fields": {
                        "idot": (0, 22), "gnss_data_info": (22, 41), "gnss_week": (41, 60), "gnss_l2p_flag": (60, 79)
                    },
                },
                7: {
                    "parser": self._parse_obs_float,
                    "fields": {
                        "sv_accuracy": (0, 22),
                        "sv_health": (22, 41),
                        "gnss_tgd_bgd": (41, 60),
                        "gnss_iodc_groupdelay": (60, 79),
                    },
                },
                8: {
                    "parser": self._parse_obs_float,
                    "fields": {
                        "transmission_time": (0, 22),
                        "gnss_interval": (22, 41),
                        # 'spare1':               (41, 60), # Blank field
                        # 'spare2':               (60, 79), # Blank field
                    },
                },
            },
        )

        return itertools.chain([header_parser], itertools.repeat(data_parser))

    def _parse_float(self, line, _):
        """Parse float entries of RINEX header to instance variable 'meta'
        """
        for k, v in line.items():
            self.meta.setdefault(self.vars["date"], dict()).update({k: _float(v)})

    def _parse_ion_alpha(self, line, _):
        """Parse entries of RINEX header `ION ALPHA` to instance variable `meta`.

                self.meta[<date>]['iono_para'] = { <gnss_id>: { 'para': [<parameter list>] }}

           Note: The ionospheric correction parameters are GNSS dependent, but we assume that is only used for GPS in
                 RINEX format version 2.11.

            ===========  =================  ======================================================================
             GNSS ID      Parameters         Description
            ===========  =================  ======================================================================
             GPSA         alpha0 - alpha3    Parameters needed for Klobuchar model (see Section 20.3.3.5.2.5 in
                                             :cite:`is-gps-200h`)
            ===========  =================  ======================================================================
        """
        para = list()
        date = self.vars["date"]

        # Save ionospheric correction parameters in a list
        for field in sorted([f for f in line if f.startswith("ion_a")]):
            para.append(_float(line[field]))

        self.meta.setdefault(date, dict())
        self.meta[date].setdefault("iono_para", dict()).update({"GPSA": {"para": para}})

    def _parse_ion_beta(self, line, _):
        """Parse entries of RINEX header `ION BETA` to instance variable `meta`.

                self.meta[<date>]['iono_para'] = { <gnss_id>: { 'para': [<parameter list>] }}

           Note: The ionospheric correction parameters are GNSS dependent, but we assume that is only used for GPS in
                 RINEX format version 2.11.

            ===========  =================  ======================================================================
             GNSS ID    Parameters               Description
            ===========  =================  ======================================================================
             GPSB         beta0 - beta3     Parameters needed for Klobuchar model (see Section 20.3.3.5.2.5 in
                                            :cite:`is-gps-200h`)
            ===========  =================  ======================================================================
        """
        para = list()
        date = self.vars["date"]

        # Save ionospheric correction parameters in a list
        for field in sorted([f for f in line if f.startswith("ion_b")]):
            para.append(_float(line[field]))

        self.meta.setdefault(date, dict())
        self.meta[date].setdefault("iono_para", dict()).update({"GPSB": {"para": para}})

    def _parse_leap_seconds(self, line, _):
        """Parse entries of RINEX header `LEAP SECONDS` to instance variable `meta`.

            self.meta[<date>]['leap_seconds'] = { 'leap_seconds': <value>,
                                                  'future_past_leap_seconds': <value>,
                                                  'week': <value>,
                                                  'week_day': <value>,
                                                  'time_sys': <system> }
        """
        # NOTE: Adding of following fields are done to be consistent with data structure of RINEX 3.xx navigation format
        #      parser.
        line.update({"future_past_leap_seconds": "", "week": "", "week_day": "", "time_sys": ""})

        date = self.vars["date"]
        self.meta.setdefault(date, dict())
        for field in line:
            self.meta[date].setdefault("leap_seconds", dict()).update({field: line[field]})

    def _parse_obs_float(self, line, cache):
        """Parse float entries of RINEX navigation data block to instance variable 'data'.
        """

        # TODO: RINEX header in between the navigation message blocks is not handled so far!!! Following lines are only
        #      a workaround to skip these RINEX header lines.
        if "skip_additional_header_line" in cache:
            return

        for k, v in line.items():
            self.data.setdefault(k, list()).append(_float(v))

    def _parse_string(self, line, _):
        """Parse string entries of RINEX header to instance variable 'meta'
        """
        for k, v in line.items():
            self.meta.setdefault(self.vars["date"], dict()).update({k: v})

    def _parse_string_list(self, line, _):
        """Parse string entries of RINEX header to instance variable 'meta' in a list
        """
        date = self.vars["date"]
        for k, v in line.items():
            self.meta.setdefault(date, dict())
            self.meta[date].setdefault(k, list()).append(v)

    def _parse_observation_epoch(self, line, cache):
        """Parse observation epoch information of RINEX navigation data record

        Only the last 2-digits of the year is given in the observation epoch, therefore it is necessary to get the
        complete 4-digit year. The year can only be determined corrrectly in the time period 1980 until 2079.
        """

        # TODO: RINEX header in between the navigation message blocks is not handled so far!!! Following lines are only
        #      a workaround to skip these RINEX header lines.
        if line["sat_clock_drift_rate"][0].isalpha():
            cache["skip_additional_header_line"] = True
            return

        # Get correct 4-digit year (in observation epoch only 2-digit year is given)
        if (int(line["year"]) >= 80) and (int(line["year"]) <= 99):
            year = int("19" + line["year"].zfill(2))
        else:
            year = int("20" + line["year"].zfill(2))

        # Create Time object
        time = (
            "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}"
            "".format(
                year=year,
                month=int(line["month"]),
                day=int(line["day"]),
                hour=int(line["hour"]),
                minute=int(line["minute"]),
                second=float(line["second"]),
            )
        )

        # Fill temporary Dataset
        self.data.setdefault("time", list()).append(time)
        sat = self.system + line["sat"].zfill(2)
        self.data.setdefault("system", list()).append(sat[0])
        self.data.setdefault("satellite", list()).append(sat)
        self.data.setdefault("sat_clock_bias", list()).append(_float(line["sat_clock_bias"]))
        self.data.setdefault("sat_clock_drift", list()).append(_float(line["sat_clock_drift"]))
        self.data.setdefault("sat_clock_drift_rate", list()).append(_float(line["sat_clock_drift_rate"]))

    def _parse_time_system_corr(self, line, _):
        """Parse entries of RINEX header `DELTA-UTC: A0,A1,T,W` to instance variable `meta`.

                self.meta[<date>]['time_sys_corr'] = { <corr_type>: { 'a0': <value>,
                                                                      'a1': <value>,
                                                                      't':  <value>,
                                                                      'w':  <value> }}

           The 'DELTA-UTC: A0,A1,T,W' header entry seems to be only defined for GPS in RINEX navigation file format
           2.11 (see :cite:`rinex2`) with following correction type:

            ===========  =====================  ======================
             Type         Parameters             Description
            ===========  =====================  ======================
             GPUT         a0, a1                 GPS to UTC
            ===========  =====================  ======================

        """
        date = self.vars["date"]
        self.meta.setdefault(date, dict())
        self.meta[date].setdefault("time_sys_corr", dict()).update(
            {"GPUT": {"a0": _float(line["a0"]), "a1": _float(line["a1"]), "t": int(line["t"]), "w": int(line["w"])}}
        )

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return [  # TODO self._check_nav_message,
            self._rename_fields_based_on_system, self._time_system_correction, self._determine_message_type
        ]

    def _rename_fields_based_on_system(self):
        """Rename general GNSS fields to GNSS specific ones

        Several fields in the RINEX navigation message have a different meaning depending on the GNSS. As a first step
        the RINEX navigation parser reads these fields in a general field variable like 'gnss_data_info',
        'gnss_interval', 'gnss_iodc_groupdelay', 'gnss_l2p_flag' and 'gnss_tgd_bgd'. As a second step the general
        fields are renamed in this routine depending on the GNSS. In the following table the relationship between
        the general field and the new GNSS dependent fields are shown:

        ===================== ==================== =================================================================
         General field         New field            Description
        ===================== ==================== =================================================================
        gnss_data_info                              Depending on GNSS this field has different meaning:
                               E: data_source        - Galileo: Data source information about the broadcast
                                                       ephemeris block, that means if the ephemeris block is based
                                                       on FNAV or INAV navigation message.
                               G: codes_l2           - GPS: Codes on L2 channel. Indication which codes are used
                                                       on L2 channel (P code, C/A code). See section 20.3.3.3.1.2
                                                       in :cite:`is-gps-200h`).
                               J: codes_l2           - QZSS: Codes on L2 channel. Indication if either C/A- or P-
                                                       code is used on L2 channel (0: spare, 1: P-code, 2: L1C/A
                                                       code). See section 4.1.2.7 in :cite:`is-qzss-pnt-001`.
                               C, I, R: None         - BeiDou, IRNSS and GLONASS: not used

        gnss_interval                               Interval indicates either the curve-fit interval for GPS or QZSS
                                                    ephemeris or for BeiDou the extrapolation interval for clock
                                                    correction parameters:
                                                    BeiDou to the clock correction parameters:
                               C: age_of_clock_corr   - BeiDou: Age of data, clock (AODC) is the extrapolated interval
                                                        of clock correction parameters. It indicates the time
                                                        difference between the reference epoch of clock correction
                                                        parameters and the last observation epoch for extrapolating
                                                        clock correction parameters. Meaning of AODC
                                                         < 25  Age of the satellite clock correction parameters in hours
                                                           25  Age of the satellite clock correction parameters is two days
                                                           ...
                                                        See section 5.2.4.9 in :cite:`bds-sis-icd`.
                               G: fit_interval        - GPS: Indicates the curve-fit interval used by the GPS Control
                                                        Segment in determining the ephemeris parameters, which is
                                                        given in HOURS (see section 6.6 in :cite:`rinex2`).
                               J: fit_interval        - QZSS: Fit interval is given as flag (see section 4.1.2.7 in
                                                        :cite:`is-qzss-pnt-001`)
                                                              0 - 2 hours
                                                              1 - more than 2 hours
                                                              blank - not known
                               E, I, R: None          - Galileo, IRNSS and GLONASS: not used

        gnss_iodc_groupdelay                        Clock issue of data or group delay depending on GNSS:
                               C: tgd_b2_b3           - BeiDou: B2/B3 TGD2
                               E: bgd_e1_e5b          - Galileo: E1-E5b BGD (see section 5.1.5 in :cite:`os-sis-icd`)
                               G: iodc                - GPS: IODC (Clock issue of data indicates changes
                                                        (set equal to IODE))
                               J: iodc                - QZSS: IODC
                               I, R: None             - IRNSS and GLONASS: not used

        gnss_l2p_flag                               L2 P-code data flag is only used by GPS and QZSS:
                               G: l2p_flag            - GPS: When bit 1 of word four is a "1", it shall indicate that
                                                        the NAV data stream was commanded OFF on the P-code of the L2
                                                        channel (see section 20.3.3.3.1.6 in :cite:`is-gps-200h`).
                               J: l2p_flag            - QZSS: L2P data flag set to 1 since QZSS does not track L2P.
                                                        See section 4.1.2.7 in :cite:`is-qzss-pnt-001`.
                               C, E, I, R: None       - BeiDou, Galileo, IRNSS and GLONASS: not used

        gnss_tgd_bgd                                Total group delay (TGD) or broadcast group delay (BGD) for
                                                    Galileo:
                               C: tgd_b1_b3           - BeiDou: B1/B3 TGD1
                               E: bgd_e1_e5a          - Galileo: E1-E5a BGD (see section 5.1.5 in :cite:`os-sis-icd`)
                               G: tgd                 - GPS: TGD (:math:`L_1 - L_2` delay correction term. See
                                                        section 20.3.3.3.3.2 in :cite:`is-gps-200h`.)
                               I: tgd                 - IRNSS: TGD
                               J: tgd                 - QZSS: TGD
                               R: None                - GLONASS: not used
        """
        for field, names in SYSNAMES.items():
            if field in self.data:
                try:
                    new_field = names[self.system]
                except KeyError:
                    # Remove blank RINEX navigation fields, which does not support this field
                    del self.data[field]
                    continue

                if field != new_field:
                    self.data[new_field] = self.data[field]
                    del self.data[field]

    def _check_nav_message(self):
        """Check correctness of navigation message.

        Issue of ephemeris data (IODE) and clock issue of data (IODC) should be equal (see Sections 20.3.3.3.1.5,
        20.3.3.4.1 and 20.3.4.4 in :cite:`is-gps-200h`). If this is not the case, a data set cutover has occurred and
        new data must be collected.
        """
        # TODO: Check has to be done GNSS dependent!!!
        # TODO: What should be done, if this happens? Use of another navigation message?
        if set(self.data["gnss_iodc_groupdelay"]).difference(set(self.data["iode"])):
            log.fatal("Issue of ephemeris data (IODE) and clock issue of data (IODC) are not equal.")

    def _determine_message_type(self):
        """Determine navigation message type and save it in 'data' dictionary under 'nav_type' key

        The navigation message type is dependent on the GNSS:
            GPS:      GPS provides LNAV and the newer CNAV and MNAV messages. LNAV and CNAV messages are civil
                      navigation messages, whereas MNAV is a military message.
            Galileo:  Galileo provides the Freely accessible (F/NAV), the Integriy (I/NAV), Commercial (C/NAV) and
                      Governmental (G/NAV) Navigation Message.

        Where can handle at the moment LNAV, F/NAV and I/NAV navigation messages. For further information see also
        in Section 8.3 in :cite:`rinex3` and section 5.1.3 in :cite:`os-sis-icd` (for Galileo F/NAV and I/NAV
        message).

        The mask for checking I/NAV and F/NAV navigation messages uses bit 8 and 9 of the data source entry in 
        RINEX 'BROADCAST ORBIT - 5' record. If bit 8 is set (e.g. 258), than the satellite clock correction is given 
        for E5a-E1 (I/NAV message) and if bit 9 (e.g. 513, 516, 517) is set, than the satellite clock correction is
        given for E5b-E1 (F/NAV message).        
        """
        nav_mask = {0b1000000000: "INAV", 0b0100000000: "FNAV"}  # Mask definition for Galileo navigation messages
        self.data["nav_type"] = np.full(len(self.data["time"]), "", dtype=object)

        for sys in set(self.data["system"]):
            idx = np.array(self.data["system"]) == sys

            if sys == "G":
                # NOTE: In Steigenberger et al. (2015): 'Performance Evaluation of the Early CNAV Navigation Message' do
                #      they use for CNAV navigation messages the RINEX 'IODE' record for the change rate of the semi-
                #      major axis 'aDot'. 'IODE' is an integer and 'aDot' a float value.
                if not np.all(_isinteger(np.array(self.data["iode"])[idx])):
                    log.fatal(
                        "GPS navigation message does not seem to be a LNAV message. Note: Where do not handle "
                        "CNAV message so far."
                    )
                self.data["nav_type"][idx] = "LNAV"

            elif sys == "E":
                # NOTE: The RINEX data source record provides information about the Galileo navigation message.
                data_source = np.array(self.data["data_source"])[idx]
                data_source[data_source == None] = 0  # Bitwise operations does not handle NoneType values
                type_tmp = np.full(len(data_source), None)

                # Loop over different Galileo navigation masks for navigation message detection
                for mask, nav_type in nav_mask.items():
                    mask_idx = np.bitwise_and(data_source.astype(int), mask).astype(bool)
                    type_tmp[mask_idx] = nav_type

                self.data["nav_type"][idx] = type_tmp

    def _time_system_correction(self):
        """Apply correction to given time system for getting GPS time scale

        Following relationship are given between GNSS time scale (either BeiDou, Galileo, IRNSS or QZSS)
        :math:`t_{GNSS}` and GPS time scale :math:`t_{GPS}` (see Section 2.1.4 in :cite:`teunissen2017`):
        .. math::
              t_{GPS}  = t_{GNSS} + \Delta t

        The time offset :math:`\Delta t` is 0 s for Galileo, IRNSS and QZSS and for BeiDou 14 s. All these time scales
        are related to the International Atomic Time (TAI) by a certain time offset. In addition the GNSS week number
        is different depending on the GNSS. Galileo, IRNSS and QZSS are referring to the same GPS week, whereas BeiDou
        week starts at GPS week 1356.

        In this routine the given navigation epoch (time of clock (toc)), time of ephemeris (toe) and the broadcast
        ephemeris transmission time will be transformed to GPS time scale for BeiDou, Galileo, IRNSS and QZSS.

        It can happen that the transmission time is referring to another GPS week as the navigation epoch (time of clock
        (toc)). The transmission time will be adjusted to the same GPS week of the navigation epoch.

        TODO: Is it necessary to adjust also toe, if given in RINEX format?
        TODO: Should this function be an editor instead?
        """
        system = self.system
        valid_systems = ["C", "E", "G", "I", "J"]  # R-GLONASS and S-SBAS is not handled so far in Where

        if system not in valid_systems:
            log.fatal(
                "Time system '{}' in file {} is not handled in Where. Time systems for following GNSS identifiers"
                " are defined: {}.",
                system,
                self.file_path,
                (", ").join(valid_systems),
            )

        # Time conversion is only necessary for C-BeiDou navigation files (Note: M-Mixed are not defined for RINEX 2.11
        # format)
        if system is "C":

            # Convert observation time entries to GPS time scale by adding system time offset
            self.data["time"] = [
                dateutil.parser.parse(t) + timedelta(seconds=SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get(system, 0))
                for t in self.data["time"]
            ]

            self.data["toe"] = [
                t + timedelta(seconds=SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get(system, 0)) for t in self.data["toe"]
            ]

            self.data["transmission_time"] = [
                t + timedelta(seconds=SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get(system, 0))
                for t in self.data["transmission_time"]
            ]

            self.data["gnss_week"] = [
                w + SYSTEM_TIME_OFFSET_TO_GPS_WEEK.get(system, 0) for w in self.data["gnss_week"]
            ]

        # Convert time data entries to Time object
        self.data["time"] = Time(self.data["time"], scale="gps")

        # TODO hjegei: Workaround -> better would if Time object can handle gpsweek as input format!!!!
        # TODO: Different GNSS week definition has to be handled here!!!!
        # self.data['toe'] = Time(self.data['gnss_week'], format='gpsweek') + TimeDelta(self.data['toe'], format='sec')
        # self.data['transmission_time'] = Time(self.data['gnss_week'], format='gpsweek') +
        #                                      TimeDelta(self.data['transmission_time'], format='sec')
        jd_day, jd_frac = gnss.gpssec2jd(np.array(self.data["gnss_week"]), np.array(self.data["toe"]))
        self.data["toe"] = Time(val=jd_day, val2=jd_frac, format="jd", scale="gps")

        jd_day, jd_frac = gnss.gpssec2jd(np.array(self.data["gnss_week"]), np.array(self.data["transmission_time"]))
        self.data["transmission_time"] = Time(val=jd_day, val2=jd_frac, format="jd", scale="gps")

        # Handling of week crossovers - refer time of ephemeris (toe) and transmission time to same GPS week as
        # navigation epoch (time of clock (toc))
        # TODO: It is only a workaround. This should be done before generating a TimeObject!!!! Use 28.02.2016 as check.
        # TODO: Is it necessary for toe? Or is toe always refered to current GPS week?
        for field in ["toe", "transmission_time"]:
            gpssec = self.data[field].gpssec.copy()
            time_diff = self.data["time"].gpssec - gpssec
            if np.any(time_diff > 302400):
                idx = time_diff > 302400
                gpssec[idx] += 604800
            elif np.any(time_diff < -302400):
                idx = time_diff < -302400
                gpssec[idx] -= 604800
            jd_day, jd_frac = gnss.gpssec2jd(np.array(self.data["gnss_week"]), gpssec)
            self.data[field] = Time(val=jd_day, val2=jd_frac, format="jd", scale="gps")

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Store GNSS RINEX navigation data in a dataset

        Args:
            dset (Dataset): The Dataset where broadcast ephemeris are stored with following fields:

            ==================== ====== ================ ===============================================================
            Field                System Unit             Description
            ==================== ====== ================ ===============================================================
            age_of_clock_corr    C                       BeiDou: Age of data, clock (AODC) is the extrapolated interval
                                                           of clock correction parameters. It indicates the time
                                                           difference between the reference epoch of clock correction
                                                           parameters and the last observation epoch for extrapolating
                                                           clock correction parameters. Meaning of AODC
                                                           < 25  Age of the satellite clock correction parameters in
                                                                 hours
                                                             25  Age of the satellite clock correction parameters is
                                                                 two days
                                                            ...
                                                           See section 5.2.4.9 in :cite:`bds-sis-icd`.
            bgd_e1_e5b            E                      Galileo: group delay E1-E5b BGD (see section 5.1.5 in
                                                           :cite:`os-sis-icd`)
            cic, cis             CEGIJ  rad              Correction coefficients of inclination
            crc, crs             CEGIJ  m                Correction coefficients of orbit radius
            cuc, cus             CEGIJ  rad              Correction coefficients of argument of perigee
            codes_l2               G J                   GPS: Codes on L2 channel. Indication which codes are used
                                                           on L2 channel (P code, C/A code). See section 20.3.3.3.1.2
                                                           in :cite:`is-gps-200h`).
                                                         QZSS: Codes on L2 channel. Indication if either C/A- or P-
                                                           code is used on L2 channel (0: spare, 1: P-code, 2: L1C/A
                                                           code). See section 4.1.2.7 in :cite:`is-qzss-pnt-001`.
            data_source           E                      Galileo: Data source information about the broadcast
                                                           ephemeris block, that means if the ephemeris block is based
                                                           on FNAV or INAV navigation message.
            delta_n              CEGIJ  rad/s            Mean motion difference from computed value
            e                    CEGIJ                   Eccentricity of the orbit
            fit_interval           G J                   GPS: Indicates the curve-fit interval used by the GPS Control
                                                          Segment in determining the ephemeris parameters, which is
                                                          given in HOURS (see section 6.6 in :cite:`rinex2`).
                                                        QZSS: Fit interval is given as flag (see section 4.1.2.7 in
                                                         :cite:`is-qzss-pnt-001`):
                                                            0 - 2 hours
                                                            1 - more than 2 hours
                                                            blank - not known
            gnss_week            CEGIJ                   Week number of ephemeris reference epoch, whichs depends on the
                                                         used GNSS. The week number is converted to GPS week.
            i0                   CEGIJ  rad              Inclination angle at the reference time
            idot                 CEGIJ  rad/s            Rate of change of inclination angle
            iodc                   G J                   GPS: IODC (Clock issue of data indicates changes (set equal to
                                                           IODE))
                                                         QZSS: IODC
            iode                 CEGIJ                   Ephemeris issue of data indicates changes to the broadcast
                                                         ephemeris:
                                                           - GPS: Ephemeris issue of data (IODE), which is set equal to
                                                             IODC
                                                           - Galileo: Issue of Data of the NAV batch (IODnav)
                                                           - QZSS: Ephemeris issue of data (IODE)
                                                           - BeiDou: Age of Data Ephemeris (AODE)
                                                           - IRNSS: Issue of Data, Ephemeris and Clock (IODEC)
            lp2_flag               G J                   L2 P-code data flag:
                                                           - GPS: When bit 1 of word four is a "1", it shall indicate
                                                             that the NAV data stream was commanded OFF on the P-code
                                                             of the L2 channel (see section 20.3.3.3.1.6 in
                                                             :cite:`is-gps-200h`).
                                                           - QZSS: L2P data flag set to 1 since QZSS does not track L2P.
                                                             See section 4.1.2.7 in :cite:`is-qzss-pnt-001`.
            m0                   CEGIJ  rad              Mean anomaly at reference epoch
            omega                CEGIJ  rad              Argument of perigee
            Omega                CEGIJ  rad              Longitude of ascending node of orbit plane at weekly epoch
            Omega_dot            CEGIJ  rad/s            Rate of change of right ascension of the ascending node
            sat_clock_bias       CEGIJ  s                Satellite clock offset from GPS time
            sat_clock_drift      CEGIJ  s/s              Satellite clock frequency offset
            sat_clock_drift_rate CEGIJ  :math:`s/s^2`    Satellite clock frequency drift
            satellite            CEGIJ                   Satellite PRN number
            sqrt_a               CEGIJ  :math:`\sqrt{m}` Square root of semi-major axis of the orbit
            sv_accuracy          CEGIJ  m                Satellite accuracy index, which is different for GNSS:
                                                           - GPS: SV accuracy in meters (see section 20.3.3.3.1.3 in
                                                             :cite:`is-gps-200h`)
                                                           - Galileo: SISA signal in space accuracy in meters
                                                             (see section 5.1.11 in :cite:`os-sis-icd`)
                                                           - BeiDou: Is that the user range accuracy index (URAI) in
                                                             meters (see section 5.2.4.5 in :cite:`bds-sis-icd`)
                                                           - QZSS: Is that the user range accuracy index?
                                                             (see table 4.1.2-4 in :cite:`is-qzss-pnt-001`)
                                                           - IRNSS: User range accuracy in meters (see section 6.2.1.4
                                                             :cite:`irnss-icd-sps`)
            sv_health            CEGIJ                   The definition of the satellite vehicle health flags depends on
                                                         GNSS:
                                                           - GPS: see section 20.3.3.3.1.4 in :cite:`is-gps-200h`
                                                           - Galileo: see section 5.1.9.3 in :cite:`os-sis-icd`
                                                           - BeiDou: see section 5.2.4.6 in :cite:`bds-sis-icd`
                                                           - QZSS: see section 4.1.2.3 in :cite:`is-qzss-pnt-001`
                                                           - IRNSS: see section 6.2.1.6 in :cite:`irnss-icd-sps`
            system               CEGIJ                   GNSS identifier
            tgd                    GIJ                   Total group delay (TGD) for GPS, IRNSS and QZSS:
                                                           - GPS: TGD (:math:`L_1 - L_2` delay correction term. See
                                                             section 20.3.3.3.3.2 in :cite:`is-gps-200h`.)
            tgd_b1_b3            C                       BeiDou: total group delay (TGD1) for frequencies B1/B3
            tgd_b2_b3            C                       BeiDou: total group delay (TGD2) for frequencies B2/B3
            time                                         Time of clock (Toc), which is related to GPS time scale. That
                                                         means all the different GNSS time systems (GPS: GPS time,
                                                         Galileo: GAL time, QZSS: QZS time, BeiDou: BDT time, IRNSS:
                                                         IRN time) are converted to GPS time scale.
            time.data                                    TODO
            time.utc                                     TODO
            toe                  CEGIJ  s                Time of ephemeris, that means fractional part of current GPS
                                                         week of ephemeris reference epoch. The week is dependent
                                                         on GNSS (GPS: GPS week, Galileo: GAL week, QZSS: GPS week,
                                                         BeiDou: BDT week, IRNSS: IRN week), therefore the different
                                                         GNSS weeks are converted to GPS week.
            transmission_time    CEGIJ  s                Transmission time of message converted to GPS time scale.
            ==================== ====== ================ ===============================================================

            and following Dataset `meta` data, whereby the `meta` data are saved in a dictionary with current days
            as keys:

            ==================== ======= ===========================================================================
             Entry                Type    Description
            ==================== ======= ===========================================================================
            comment               list    List with comment lines
            file_created          str     Date of file creation
            file_type             str     File type
            iono_para             dict    Dictionary with GNSS dependent ionospheric correction parameters (assumed
                                          that is defined only for GPS by using RINEX 2.11 format)
            leap_seconds          dict    Dictionary with information related to leap seconds (compatible with RINEX
                                          3.xx navigation file data struture)
            program               str     Name of program creating current file
            run_by                str     Name of agency creating current file
            time_sys_corr         dict    Dictionary with GNSS time system corrections
            version               str     Format version
            ==================== ======= ===========================================================================
        """

        dset.num_obs = len(self.data["time"])
        dset.meta.update(self.meta)

        for k, v in self.data.items():
            if k in ["time", "toe", "transmission_time"]:
                dset.add_time(k, val=v, scale="gps")
            elif k in ["nav_type", "satellite", "system"]:
                dset.add_text(k, val=v)
            else:
                if isinstance(v, list):
                    dset.add_float(k, val=np.array(v))


# TODO: Maybe better to have this routine in a module.
def _float(value):
    """Convert string to float value

    Convert a string to a floating point number (including, e.g. -0.5960D-01). Whitespace or empty value is set to 0.0.

    Args:
        value (str):   string value

    Returns:
        float: float value
    """
    if value.isspace() or not value:
        return 0.0
    else:
        return float(value.replace("D", "e"))


# TODO: Maybe better to have this routine in a module.
def _isinteger(x):
    """Check for integer type in given array

    Args:
        x (float or numpy.ndarray):  Float numbers

    Returns:
        numpy.ndarray:  Array boolean values, whereby True means integer value
    """
    return np.equal(np.mod(x, 1), 0)
