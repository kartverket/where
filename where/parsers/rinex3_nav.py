"""A parser for reading GNSS RINEX v3.03 navigation file (exception GLONASS and SBAS)

Example:
--------

    from where import data
    from where import parsers

    # Parse data
    parser = parsers.parse('rinex3_nav', rundate=rundate)

    # Create a empty Dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=0, empty=True)

    # Fill Dataset with parsed data
    parser.write_to_dataset(dset)

Description:
------------

Reads GNSS data from files in the RINEX navigation file format 3.03 (see :cite:`rinex3`). Note that beside the given
day also the navigation files from the day before and the day after is read. An exception is also, that this parser
does not handle GLONASS and SBAS navigation messages. All navigation time epochs (time of clock (toc)) are
converted to GPS time scale.


"""

# Standard library imports
from datetime import timedelta
import dateutil.parser
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers import parser
from where.lib import config
from where.lib import files
from where.lib import gnss
from where.lib import log
from where.lib.time import Time, TimeDelta

# TODO: SYSTEM_TIME_OFFSET_TO_GPS_SECOND & SYSTEM_TIME_OFFSET_TO_GPS_WEEK should be placed in constans.conf

# The constant shows the 'second' relation between the GNSS time systems to the GPS time scale. Galileo (E),
# QZSS (I) and IRNSS (J) uses the same second as the GPS time scale, but BeiDou (C) BDT time system is 14 seconds behind
# GPS time.
SYSTEM_TIME_OFFSET_TO_GPS_SECOND = dict(C=14, E=0, G=0, I=0, J=0)

# The constant shows the relation between the GNSS week to the GPS week. Galileo (E), IRNSS (I) and QZSS (J) week
# corresponds to GPS week, whereas BeiDou (C) week starts at GPS week 1356.
SYSTEM_TIME_OFFSET_TO_GPS_WEEK = dict(C=1356, E=0, G=0, I=0, J=0)

SYSNAMES = dict(
    gnss_data_info={"G": "codes_l2", "J": "codes_l2", "E": "data_source"},
    gnss_interval={"G": "fit_interval", "J": "fit_interval", "C": "age_of_clock_corr"},
    gnss_iodc_groupdelay={"G": "iodc", "J": "iodc", "E": "bgd_e1_e5b", "C": "tgd_b2_b3"},
    gnss_l2p_flag={"G": "l2p_flag", "J": "l2p_flag"},
    gnss_tgd_bgd={"G": "tgd", "J": "tgd", "E": "bgd_e1_e5a", "C": "tgd_b1_b3", "I": "tgd"},
)


@plugins.register
class Rinex3NavParser(parser.Parser):
    """A parser for reading RINEX navigation file

    The parser reads GNSS broadcast ephemeris in RINEX format 3.03 (see :cite:`rinex3`) except for GLONASS and SBAS.

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
        _parse_file()             Read a data file and parse the content
        _parse_float()            Parse float entries of SP3 header to instance variable 'meta'
        _parse_ionospheric_corr() Parse entries of RINEX header `IONOSPHERIC CORR` to instance variable `meta`.
        _parse_leap_seconds()     Parse entries of RINEX header `LEAP SECONDS` to instance variable `meta`.
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

        Overridden to read data for three days, the given date and the day before.

        TODO: Is that necessary to read three consecutive days of broadcast ephemeris???? No interpolation needed. The
              valid RINEX epoch is used for given time. The routine get_brdc_idx(dset) has to be changed.
        """
        date_to_read = self.rundate - timedelta(days=1)

        while date_to_read <= (self.rundate + timedelta(days=1)):
            self.vars.update(config.date_vars(date_to_read))
            self.vars["date"] = "{yyyy}-{mm}-{dd}".format(**self.vars)
            file_path = files.path(self.file_key, file_vars=self.vars)
            self.dependencies.append(file_path)

            with files.open_path(file_path, mode="rt") as fid:
                self.parse_file(fid)

            date_to_read += timedelta(days=1)

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
                #      3.03           N: GNSS NAV DATA    E: GALILEO          RINEX VERSION / TYPE
                "RINEX VERSION / TYPE": {
                    "parser": self._parse_string,
                    "fields": {"version": (0, 20), "file_type": (20, 21), "sat_sys": (40, 41)},
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
                # BDSA    .1397E-07   .0000E+00  -.5960E-07   .5960E-07       IONOSPHERIC CORR
                # BDSB    .1106E+06  -.3277E+05  -.2621E+06   .1966E+06       IONOSPHERIC CORR
                "IONOSPHERIC CORR": {
                    "parser": self._parse_ionospheric_corr,
                    "fields": {
                        "gnss_id": (0, 4),
                        "para_1": (5, 17),
                        "para_2": (17, 29),
                        "para_3": (29, 41),
                        "para_4": (41, 53),
                        "time_mark": (54, 55),
                        "sv_id": (56, 58),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # BDUT -5.5879354477e-09-0.000000000e+00     14 1886          TIME SYSTEM CORR
                # GAUT  0.0000000000e+00 0.000000000e+00 172800 1886          TIME SYSTEM CORR
                "TIME SYSTEM CORR": {
                    "parser": self._parse_time_system_corr,
                    "fields": {"corr_type": (0, 4), "a0": (5, 22), "a1": (22, 38), "t": (38, 45), "w": (45, 50)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     16    17  1851     3                                    LEAP SECONDS
                "LEAP SECONDS": {
                    "parser": self._parse_leap_seconds,
                    "fields": {
                        "leap_seconds": (0, 6),
                        "future_past_leap_seconds": (6, 12),
                        "week": (12, 18),
                        "week_day": (18, 24),
                        "time_sys": (24, 27),
                    },
                },
            },
        )

        # Parser for RINEX navigation data blocks

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
        # E11 2016 02 28 22 00 00  .654120231047E-04  .109707798401E-10  .000000000000E+00
        #       .400000000000E+01  .129375000000E+02  .319691887879E-08 -.292934515480E+01
        #       .460073351860E-06  .329698785208E-03  .683590769768E-05  .544061308098E+04
        #       .792000000000E+05  .391155481339E-07 -.125937621871E+01  .316649675369E-07
        #       .967916388734E+00  .197406250000E+03 -.653087089047E+00 -.571166648526E-08
        #       .276797244002E-09  .257000000000E+03  .188600000000E+04  .000000000000E+00
        #      -.100000000000E+01  .000000000000E+00 -.249128788710E-07 -.225845724344E-07
        #       .798850000000E+05  .000000000000E+00  .000000000000E+00  .000000000000E+00
        data_parser = parser.define_parser(
            end_marker=lambda _l, _ln, next_line: next_line[
                0:1
            ].isalpha(),  # Observation line start with satellite system identifier
            label=lambda _l, line_num: line_num,
            parser_def={
                1: {
                    "parser": self._parse_observation_epoch,
                    "fields": {
                        "system": (0, 1),
                        "sat_num": (1, 3),
                        "year": (4, 8),
                        "month": (9, 11),
                        "day": (12, 14),
                        "hour": (15, 17),
                        "minute": (18, 20),
                        "second": (21, 23),
                        "sat_clock_bias": (23, 42),
                        "sat_clock_drift": (42, 61),
                        "sat_clock_drift_rate": (61, 80),
                    },
                },
                2: {
                    "parser": self._parse_obs_float,
                    "fields": {"iode": (0, 23), "crs": (23, 42), "delta_n": (42, 61), "m0": (61, 80)},
                },
                3: {
                    "parser": self._parse_obs_float,
                    "fields": {"cuc": (0, 23), "e": (23, 42), "cus": (42, 61), "sqrt_a": (61, 80)},
                },
                4: {
                    "parser": self._parse_obs_float,
                    "fields": {"toe": (0, 23), "cic": (23, 42), "Omega": (42, 61), "cis": (61, 80)},
                },
                5: {
                    "parser": self._parse_obs_float,
                    "fields": {"i0": (0, 23), "crc": (23, 42), "omega": (42, 61), "Omega_dot": (61, 80)},
                },
                6: {
                    "parser": self._parse_obs_float,
                    "fields": {
                        "idot": (0, 23),
                        "gnss_data_info": (23, 42),
                        "gnss_week": (42, 61),
                        "gnss_l2p_flag": (61, 80),
                    },
                },
                7: {
                    "parser": self._parse_obs_float,
                    "fields": {
                        "sv_accuracy": (0, 23),
                        "sv_health": (23, 42),
                        "gnss_tgd_bgd": (42, 61),
                        "gnss_iodc_groupdelay": (61, 80),
                    },
                },
                8: {
                    "parser": self._parse_obs_float,
                    "fields": {
                        "transmission_time": (0, 23),
                        "gnss_interval": (23, 42),
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

    def _parse_ionospheric_corr(self, line, _):
        """Parse entries of RINEX header `IONOSPHERIC CORR` to instance variable `meta`.

                self.meta[<date>]['iono_para'] = { <gnss_id>: { 'para': [<parameter list>],
                                                                'sv_id': <id>,
                                                                'time_mark': <value> }}

           The ionospheric correction parameters are GNSS dependent. In the following a list with GNSS idendifiers
           and belonging ionosphere correction parameters is shown:

            ===========  =================  ======================================================================
             GNSS ID      Parameters         Description
            ===========  =================  ======================================================================
             GAL          ai0 - ai2         Parameters needed for NeQuick model (see Section 5.1.6 in
                                            :cite:`os-sis-icd`)
             GPSA         alpha0 - alpha3   Parameters needed for Klobuchar model (see Section 20.3.3.5.2.5 in
                                            :cite:`is-gps-200h`)
             GPSB         beta0 - beta3     Parameters needed for Klobuchar model (see Section 20.3.3.5.2.5 in
                                            :cite:`is-gps-200h`)
             QZSA         alpha0 - alpha3   Parameters needed for Klobuchar model (see Section 5.2.4.7 in
                                            :cite:`bds-sis-icd-20`)
             QZSB         beta0 - beta3     Parameters needed for Klobuchar model (see :cite:`is-qzss`)
             BDSA         alpha0 - alpha3   Parameters needed for Klobuchar model (see :cite:`is-qzss`)
             BDSB         beta0 - beta3     Parameters needed for Klobuchar model (see Section 5.2.4.7 in
                                            :cite:`bds-sis-icd-20`)
             IRNA         alpha0 - alpha3   Parameters needed for Klobuchar model (see Appendix H in
                                            :cite:`irnss-icd-sps`)
             IRNB         beta0 - beta3     Parameters needed for Klobuchar model (see Appendix H in
                                            :cite:`irnss-icd-sps`)
            ===========  =================  ======================================================================

        """
        para = list()
        date = self.vars["date"]

        # Save ionospheric correction parameters in a list
        for field in sorted([f for f in line if f.startswith("para_")]):
            para.append(_float(line[field]))

        self.meta.setdefault(date, dict())
        self.meta[date].setdefault("iono_para", dict()).update(
            {line["gnss_id"]: {"para": para, "sv_id": line["sv_id"], "time_mark": line["time_mark"]}}
        )

    def _parse_leap_seconds(self, line, _):
        """Parse entries of RINEX header `LEAP SECONDS` to instance variable `meta`.

            self.meta[<date>]['leap_seconds'] = { 'leap_seconds': <value>,
                                                  'future_past_leap_seconds': <value>,
                                                  'week': <value>,
                                                  'week_day': <value>,
                                                  'time_sys': <system> }
        """
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

        # TODO: SBAS and GLONASS satellites has different number of input parameters. So far it is not handled in Where.
        if cache["system"] == "S" or cache["system"] == "R":
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
        """

        # TODO: RINEX header in between the navigation message blocks is not handled so far!!! Following lines are only
        #      a workaround to skip these RINEX header lines.
        if line["sat_clock_drift"][-1].isalpha():
            cache["skip_additional_header_line"] = True
            return

        # TODO: SBAS and GLONASS satellites has different number of input parameters. So far it is not handled in Where.
        cache["system"] = line["system"]
        if line["system"] == "S" or line["system"] == "R":
            return

        # Create Time object
        time = "{year}-{month:02d}-{day:02d}T{hour:02d}:{minute:02d}:{second:010.7f}" "".format(
            year=int(line["year"]),
            month=int(line["month"]),
            day=int(line["day"]),
            hour=int(line["hour"]),
            minute=int(line["minute"]),
            second=float(line["second"]),
        )

        # Fill temporary Dataset
        self.data.setdefault("time", list()).append(time)
        self.data.setdefault("system", list()).append(line["system"])
        self.data.setdefault("satellite", list()).append(line["system"] + line["sat_num"].zfill(2))
        self.data.setdefault("sat_clock_bias", list()).append(_float(line["sat_clock_bias"]))
        self.data.setdefault("sat_clock_drift", list()).append(_float(line["sat_clock_drift"]))
        self.data.setdefault("sat_clock_drift_rate", list()).append(_float(line["sat_clock_drift_rate"]))

    def _parse_time_system_corr(self, line, _):
        """Parse entries of RINEX header `TIME SYSTEM CORR` to instance variable `meta`.

                self.meta[<date>]['time_sys_corr'] = { <corr_type>: { 'a0': <value>,
                                                                      'a1': <value>,
                                                                      't':  <value>,
                                                                      'w':  <value> }}

           The time system correction parameters are dependent on the given correction type, which is shown below:

            ===========  =====================  ======================================================================
             Type         Parameters             Description
            ===========  =====================  ======================================================================
             BDUT         a0=A_0UTC, a1=A_1UTC   BDS to UTC
             GAUT         a0, a1                 GAL to UTC
             GLGP         a0=-TauC, a1=0         GLO to GPS
             GLUT         a0=-TauC, a1=0         GLO to UTC
             GPGA         a0=A0G, a1=A1G         GPS to GAL
             GPUT         a0, a1                 GPS to UTC
             IRGP         a0=A_0UTC, a1=A_1UTC   IRN to GPS
             IRUT         a0=A_0UTC, a1=A_1UTC   IRN to UTC
             QZGP         a0, a1                 QZS to GPS
             QZUT         a0, a1                 QZS to UTC
             SBUT         a0, a1                 SBAS to UTC
            ===========  =====================  ======================================================================

        """
        date = self.vars["date"]
        self.meta.setdefault(date, dict())
        self.meta[date].setdefault("time_sys_corr", dict()).update(
            {
                line["corr_type"]: {
                    "a0": _float(line["a0"]),
                    "a1": _float(line["a1"]),
                    "t": int(line["t"]),
                    "w": int(line["w"]),
                }
            }
        )

    #
    # SETUP CALCULATION
    #
    # TODO: Due to :cite:`rinex3`, was the implementation of the fit interval field from the GPS navigation message
    #       an issue for the RINEX 3.02. Some implementations wrote the flag and others wrote a time interval.
    #       The RINEX 3.03 release specifies that the fit interval should be a time period for GPS and a flag
    #       for QZSS. TPP navigation files write a flag instead of a time interval for GPS, whereby 0 = 4h and
    #       1 > 4h. Should it be handled in the RINEX parser? So far it is handled in 'remover'
    #       ignore_unavailable_orbit.py.
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return [  # TODO self.check_nav_message,
            self._rename_fields_based_on_system,
            self._time_system_correction,
            self._determine_message_type,
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
            # Invert SYSNAMES-mapping to point from system to field name
            new_fields = dict()
            for system, new_field in names.items():
                new_fields.setdefault(new_field, list()).append(system)

            # Add values to new fields based on system, and delete old field
            for new_field, systems in new_fields.items():
                self.data[new_field] = [
                    v if s in systems else None for s, v in zip(self.data["system"], self.data[field])
                ]
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
        # TODO: Does it exist a better solution to detect INAV and FNAV message types?
        # Mask definition for Galileo navigation messages
        nav_mask = {
            0b1000000001: "INAV_E1",  # = 513  <- normally used
            0b1000001001: "INAV_E1",  # = 521
            0b1000010001: "INAV_E1",  # = 529
            0b1000011001: "INAV_E1",  # = 537
            0b1000000100: "INAV_E5b",  # = 516  <- normally used
            0b1000001100: "INAV_E5b",  # = 524
            0b1000010100: "INAV_E5b",  # = 532
            0b1000011100: "INAV_E5b",  # = 540
            0b1000000101: "INAV_E1E5b",  # = 517  <- normally used
            0b1000001101: "INAV_E1E5b",  # = 525
            0b1000010101: "INAV_E1E5b",  # = 533
            0b1000011101: "INAV_E1E5b",  # = 541
            0b0100000010: "FNAV_E5a",  # = 258  <- normally used
            0b0100001010: "FNAV_E5a",  # = 266
            0b0100010010: "FNAV_E5a",  # = 274
            0b0100011010: "FNAV_E5a",  # = 282
        }

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
                    # mask_idx = np.bitwise_and(data_source.astype(int), mask).astype(bool)
                    mask_idx = data_source.astype(int) == mask
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

        TODO: Should this function be an editor instead?
        """
        system = self.meta[self.rundate.strftime("%Y-%m-%d")]["sat_sys"]

        valid_systems = ["C", "E", "G", "I", "J", "M"]  # R-GLONASS and S-SBAS is not handled so far in Where

        if system not in valid_systems:
            log.fatal(
                "Time system '{}' in file {} is not handled in Where. Time systems for following GNSS identifiers"
                " are defined: {}.",
                system,
                self.file_path,
                (", ").join(valid_systems),
            )

        # Time conversion is only necessary for M-Mixed or C-BeiDou navigation files
        if system is "M" or system is "C":
            self.data["time"] = [
                dateutil.parser.parse(t) + timedelta(seconds=SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get(s, 0))
                for s, t in zip(self.data["system"], self.data["time"])
            ]

            self.data["toe"] = [
                t + SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get(s, 0) for s, t in zip(self.data["system"], self.data["toe"])
            ]

            self.data["transmission_time"] = [
                t + SYSTEM_TIME_OFFSET_TO_GPS_SECOND.get(s, 0)
                for s, t in zip(self.data["system"], self.data["transmission_time"])
            ]

            self.data["gnss_week"] = [
                w + SYSTEM_TIME_OFFSET_TO_GPS_WEEK.get(s, 0)
                for s, w in zip(self.data["system"], self.data["gnss_week"])
            ]

        # Convert time data entries to Time object
        self.data["time"] = Time(self.data["time"], scale="gps")

        # TODO hjegei: Workaround -> better would if Time object can handle gpsweek as input format!!!!
        # self.data['toe'] = Time(self.data['gnss_week'], format='gpsweek') + TimeDelta(self.data['toe'], format='sec')
        # self.data['transmission_time'] = Time(self.data['gnss_week'], format='gpsweek') + TimeDelta(self.data['transmission_time'], format='sec')
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
                                                         See RINEX 3 format description :cite:`rinex3`.
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
            nav_type              EG                     Navigation message type (Note: At the moment only implemented
                                                         for GPS and Galileo.)
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
            iono_para             dict    Dictionary with GNSS dependent ionospheric correction parameters
            leap_seconds          dict    Dictionary with information related to leap seconds
            program               str     Name of program creating current file
            run_by                str     Name of agency creating current file
            sat_sys               str     Satellite system
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
