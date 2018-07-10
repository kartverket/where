"""A parser for reading ANTEX format 1.4 data

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='gnss_antex', file_path='igs14.atx')
    data = p.as_dict()

Description:
------------

Reads data from files in the GNSS Antenna Exchange (ANTEX) file format version 1.4 (see :cite:`antex`).


"""

# Standard library imports
import datetime
import itertools

# External library imports
import numpy as np

# Where imports
from where.parsers._parser_chain import ChainParser, ParserDef
from where.lib import log
from where.lib import plugins
from where.lib.unit import unit

# TODO: Not implemented so far in Midgard
# from midgard.config.unit import unit
# from midgard.dev import log
# from midgard.dev import plugins
# from midgard.parsers._parser_chain import ChainParser, ParserDeg


@plugins.register
class AntexParser(ChainParser):
    """A parser for reading ANTEX file

    The parser reads GNSS ANTEX format 1.4 (see :cite:`antex`).


    The 'data' attribute is a dictionary with GNSS satellite PRN or receiver antenna as key. The GNSS satellite
    antenna corrections are time dependent and saved with "valid from" datetime object entry. The dictionary looks like:

        dout = { <prn> : { <valid from>: { cospar_id:   <value>,
                                           sat_code:    <value>,
                                           sat_type:    <value>,
                                           valid_until: <value>,
                                           azimuth:     <list with azimuth values>,
                                           elevation:   <list with elevation values>,
                                           <frequency>: { azi: [<list with azimuth-elevation dependent corrections>],
                                                          neu: [north, east, up],
                                                          noazi: [<list with elevation dependent corrections>] }}},

                 <receiver antenna> : { azimuth:     <list with azimuth values>,
                                        elevation:   <list with elevation values>,
                                        <frequency>: { azi: [<array with azimuth-elevation dependent corrections>],
                                                       neu: [north, east, up],
                                                       noazi: [<list with elevation dependent corrections>] }}}

    with following entries:

    =================== =================== ========================================================================
     Value               Type                Description
    =================== =================== ========================================================================
     azi                 numpy.ndarray       Array with azimuth-elevation dependent antenna correction in [mm] with
                                             the shape: number of azimuth values x number of elevation values.
     azimuth             numpy.ndarray       List with azimuth values in [rad] corresponding to antenna corrections
                                             given in `azi`.
     cospar_id           str                 COSPAR ID <yyyy-xxxa>: yyyy -> year when the satellite was put in
                                             orbit, xxx -> sequential satellite number for that year, a -> alpha
                                             numeric sequence number within a launch
     elevation           numpy.ndarray       List with elevation values in [rad] corresponding to antenna
                                             corrections given in `azi` or `noazi`.
     <frequency>         str                 Frequency identifier (e.g. G01 - GPS L1)
     neu                 list                North, East and Up eccentricities in [m]. The eccentricities of the
                                             mean antenna phase center is given relative to the antenna reference
                                             point (ARP) for receiver antennas or to the center of mass of the
                                             satellite in X-, Y- and Z-direction.
     noazi               numpy.ndarray       List with elevation dependent (non-azimuth-dependent) antenna
                                             correction in [mm].
     <prn>               str                 Satellite code e.g. GPS PRN, GLONASS slot or Galileo SVID number
     <receiver antenna>  str                 Receiver antenna name together with radome code
     sat_code            str                 Satellite code e.g. GPS SVN, GLONASS number or Galileo GSAT number
     sat_type            str                 Satellite type (e.g. BLOCK IIA)
     valid_from          datetime.datetime   Start of validity period of satellite in GPS time
     valid_until         datetime.datetime   End of validity period of satellite in GPS time
    =================== =================== ========================================================================


    The 'meta' attribute is a dictionary with following entries:

        =================== ======= ========================================================================
         Value               Type    Description
        =================== ======= ========================================================================
         comment             list    Header commments given in list line by line
         pcv_type            str     Phase center variation type
         ref_antenna         str     Reference antenna type for relative antenna
         ref_serial_num      str     Serial number of the reference antenna
         sat_sys             str     Satellite system
         version             str     Format version
        =================== ======= ========================================================================

    Attributes:
        data (dict):                    Contains the (observation) data read from file.
        data_available (bool):          Indicator of whether data are available.
        file_path (pathlib.PosixPath):  File path.
        parser_name (str):              Parser name.
        meta (dict):                    Contains metainformation read from file.
    """
    #
    # PARSERS
    #
    def setup_parser(self):
        """Parsers defined for reading ANTEX file line by line.

           First the ANTEX header information are read and afterwards the ANTEX corrections.
        """
        # Parser for ANTEX header
        header_parser = ParserDef(
            end_marker=lambda line, _ln, _n: line[60:73] == "END OF HEADER",
            label=lambda line, _ln: line[60:].strip(),
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      1.4            M                                       ANTEX VERSION / SYST
                "ANTEX VERSION / SYST": {
                    "parser": self.parse_string, "fields": {"version": (0, 8), "sat_sys": (20, 21)}
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # A                                                           PCV TYPE / REFANT
                "PCV TYPE / REFANT": {
                    "parser": self.parse_string,
                    "fields": {"pcv_type": (0, 1), "ref_antenna": (20, 40), "ref_serial_num": (40, 60)},
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # Compiled by Arturo Villiger (AIUB),                         COMMENT
                "COMMENT": {"parser": self.parse_comment, "fields": {"comment": (0, 60)}},
            },
        )

        # Parser for ANTEX observation blocks

        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9
        #                                                             START OF ANTENNA
        # BLOCK IIA           G01                 G032      1992-079A TYPE / SERIAL NO
        #                                              0    29-JAN-17 METH / BY / # / DATE
        #      0.0                                                    DAZI
        #      0.0  17.0   1.0                                        ZEN1 / ZEN2 / DZEN
        #      2                                                      # OF FREQUENCIES
        #   1992    11    22     0     0    0.0000000                 VALID FROM
        #   2008    10    16    23    59   59.9999999                 VALID UNTIL
        # IGS14_1949                                                  SINEX CODE
        #    G01                                                      START OF FREQUENCY
        #     279.00      0.00   2319.50                              NORTH / EAST / UP
        #    NOAZI   -0.80   -0.90   -0.90   -0.80   -0.40    0.20    0.80    1.30    1.40    1.20    0.70    0.00   -0.40   -0.70   -0.90   -0.90   -0.90   -0.90
        #    G01                                                      END OF FREQUENCY
        #    G02                                                      START OF FREQUENCY
        #     279.00      0.00   2319.50                              NORTH / EAST / UP
        #    NOAZI   -0.80   -0.90   -0.90   -0.80   -0.40    0.20    0.80    1.30    1.40    1.20    0.70    0.00   -0.40   -0.70   -0.90   -0.90   -0.90   -0.90
        #    G02                                                      END OF FREQUENCY
        #                                                            END OF ANTENNA
        corr_parser = ParserDef(
            end_marker=lambda line, _ln, _: line[60:74] == "END OF ANTENNA",
            skip_line=lambda line: not line,
            label=lambda line, _ln: line[60:].strip()
            if (line[60:61].isalpha() or line[60:61] == "#")
            else "CORRECTION",
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # BLOCK IIA           G01                 G032      1992-079A TYPE / SERIAL NO
                "TYPE / SERIAL NO": {
                    "parser": self.parse_section_string,
                    "fields": {
                        "antenna_type": (0, 20), "antenna_code": (20, 40), "sat_code": (40, 50), "cospar_id": (50, 60)
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #                                              0    29-JAN-17 METH / BY / # / DATE
                # TODO 'METH / BY / # / DATE':
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      5.0                                                    DAZI
                "DAZI": {"parser": self.parse_section_float, "fields": {"dazi": (2, 8)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      0.0  17.0   1.0                                        ZEN1 / ZEN2 / DZEN
                "ZEN1 / ZEN2 / DZEN": {
                    "parser": self.parse_section_float, "fields": {"zen1": (2, 8), "zen2": (8, 14), "dzen": (14, 20)}
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #      4                                                      # OF FREQUENCIES
                "# OF FREQUENCIES": {"parser": self.parse_num_of_frequencies, "fields": {"num_freq": (0, 6)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #   1992    11    22     0     0    0.0000000                 VALID FROM
                "VALID FROM": {
                    "parser": self.parse_valid_from,
                    "fields": {
                        "year": (0, 6),
                        "month": (6, 12),
                        "day": (12, 18),
                        "hour": (18, 24),
                        "minute": (24, 30),
                        "second": (30, 43),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #   2008    10    16    23    59   59.9999999                 VALID UNTIL
                "VALID UNTIL": {
                    "parser": self.parse_valid_until,
                    "fields": {
                        "year": (0, 6),
                        "month": (6, 12),
                        "day": (12, 18),
                        "hour": (18, 24),
                        "minute": (24, 30),
                        "second": (30, 43),
                    },
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # IGS14_1949                                                  SINEX CODE
                # TODO: 'SINEX CODE':
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # ATTENTION! ROUNDED BLOCK MEAN Z-OFFSET VALUE!               COMMENT
                # TODO: 'COMMENT':
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #    G02                                                      START OF FREQUENCY
                "START OF FREQUENCY": {"parser": self.parse_section_string, "fields": {"frequency_code": (3, 6)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #     279.00      0.00   2289.30                              NORTH / EAST / UP
                "NORTH / EAST / UP": {
                    "parser": self.parse_section_float, "fields": {"north": (0, 10), "east": (10, 20), "up": (20, 30)}
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #    G02                                                      END OF FREQUENCY
                "END OF FREQUENCY": {"parser": self.save_correction, "fields": {"frequency_code": (3, 6)}},
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                #    NOAZI   +0.00   -0.17   -0.66   -1.37   -2.18   -3.03   -3.85   -4.63   -5.28   ...
                #      0.0   +0.00   -0.26   -0.82   -1.55   -2.34   -3.10   -3.82   -4.50   -5.10   ...
                #      5.0   +0.00   -0.26   -0.83   -1.56   -2.35   -3.11   -3.83   -4.50   -5.08   ...
                "CORRECTION": {
                    "parser": self.parse_correction,
                    "fields": {"values": (0, None)},  # 'None' indicates, that line is sliced until the end
                },  # of line.
            },
        )

        return itertools.chain([header_parser], itertools.repeat(corr_parser))

    #
    # HEADER PARSERS
    #
    def parse_comment(self, line, _):
        """Parse comment lines in ANTEX header.
        """
        self.meta.setdefault("comment", list()).append(line["comment"])

    def parse_string(self, line, _):
        """Parse string entries of ANTEX header.
        """
        self.parse_default_meta({k: v for k, v in line.items()}, _)

    #
    # ANTENNA CORRECTION PARSERS
    #
    def parse_section_float(self, line, cache):
        """Parse float entries of ANTEX header.
        """
        cache.update({k: float(v) for k, v in line.items()})

    def parse_section_integer(self, line, cache):
        """Parse integer entries of ANTEX header.
        """
        cache.update({k: int(v) for k, v in line.items()})

    def parse_section_string(self, line, cache):
        """Parse string entries of ANTEX header.
        """
        cache.update({k: v for k, v in line.items()})

    def parse_correction(self, line, cache):
        """Parse antenna corrections entries of ANTEX antenna section.
        """
        line = line["values"].split()

        if line[0] == "NOAZI":
            del line[0]
            cache["noazi"] = line
        else:
            del line[0]
            cache.setdefault("azi", list()).append(line)

    def parse_num_of_frequencies(self, line, cache):
        """Parse '# OF FREQUENCIES' entry of ANTEX antenna section.
        """
        cache["num_freq"] = line["num_freq"]
        cache["num_freq_counter"] = 0

    def parse_valid_from(self, line, cache):
        """Parse 'VALID FROM' entries of ANTEX antenna section.
        """
        dt = datetime.datetime(
            int(line["year"]), int(line["month"]), int(line["day"]), int(line["hour"]), int(line["minute"])
        )
        dt = dt + datetime.timedelta(float(line["second"]))
        cache["valid_from"] = dt

    def parse_valid_until(self, line, cache):
        """Parse 'VALID UNTIL' entries of ANTEX antenna section.
        """
        dt = datetime.datetime(
            int(line["year"]), int(line["month"]), int(line["day"]), int(line["hour"]), int(line["minute"])
        )
        dt = dt + datetime.timedelta(float(line["second"]))
        cache["valid_until"] = dt

    def parse_default_meta(self, line, _):
        """Add the contents of line to meta

        Args:
            line: Dict containing the fields of a line.
        """
        self.meta.update(line)

    #
    # SAVE ANTENNA CORRECTION
    #
    def save_correction(self, line, cache):
        """Save antenna correction in data structures.

        The antenna corrections are saved after reading of corrections for one frequency. Antenna correction data are
        saved in following data structure, whereby satellite antenna corrections are time dependent:

            self.data = { <prn> : { <valid from>: { cospar_id:   <value>,
                                                    sat_code:    <value>,
                                                    sat_type:    <value>,
                                                    valid_until: <value>,
                                                    azimuth:     <list with azimuth values>,
                                                    elevation:   <list with elevation values>,
                                                    <frequency>: { azi: [<list with azimuth-elevation dependent corrections>],
                                                                   neu: [north, east, up],
                                                                   noazi: [<list with elevation dependent corrections>] }}},

                          <receiver antenna> : { azimuth:     <list with azimuth values>,
                                                 elevation:   <list with elevation values>,
                                                 <frequency>: { azi: [<array with azimuth-elevation dependent corrections>],
                                                                neu: [north, east, up],
                                                                noazi: [<list with elevation dependent corrections>] }}
                        }

        """
        ant = cache["antenna_code"] if cache["sat_code"] else cache["antenna_type"]
        self.data.setdefault(ant, dict())
        freq = cache["frequency_code"]
        if "valid_from" in cache:
            dt = cache["valid_from"]
        tmp = dict()  # Temporary dictionary, where antenna correction for one frequency is saved.
        tmp[freq] = dict()

        # Save general information of ANTEX antenna section (NOTE: Has to be done only once.)
        if cache["num_freq_counter"] == 0:

            if cache["sat_code"]:  # only necessary for satellites
                if dt in self.data[ant]:
                    valid_from = "{:4d}-{:02d}-{:02d}".format(dt.year, dt.month, dt.day)
                    log.fatal(
                        f"Antenna correction for satellite PRN {ant} (SVN {cache['sat_code']}) valid from {valid_from} is not unique."
                    )

                tmp["cospar_id"] = cache["cospar_id"]
                tmp["sat_code"] = cache["sat_code"]
                tmp["sat_type"] = cache["antenna_type"]
                if "valid_until" in cache:
                    tmp["valid_until"] = cache["valid_until"]
                else:
                    tmp["valid_until"] = datetime.datetime.now()

            # Determine elevation list
            if cache["dzen"] != 0.0:
                tmp["elevation"] = np.arange(
                    90.0 - cache["zen1"], 90.0 - (cache["zen2"] + cache["dzen"]), -cache["dzen"]
                )
                tmp["elevation"] = np.radians(tmp["elevation"])

            # Determine azimuth list
            if cache["dazi"] != 0.0:
                tmp["azimuth"] = np.arange(0, 360 + cache["dazi"], cache["dazi"])
                tmp["azimuth"] = np.radians(tmp["azimuth"])

        # Save frequency dependent antenna corrections
        tmp[freq]["neu"] = [
            cache["north"] * unit.millimeter2meter,
            cache["east"] * unit.millimeter2meter,
            cache["up"] * unit.millimeter2meter,
        ]
        tmp[freq]["noazi"] = np.array(cache["noazi"])
        if "azi" in cache:
            tmp[freq]["azi"] = np.array(cache["azi"])

        # Save satellite antenna correction in data structure
        if cache["sat_code"]:
            self.data[ant].setdefault(dt, dict())

            if freq in self.data[ant][dt]:
                valid_from = "{:4d}-{:02d}-{:02d}".format(dt.year, dt.month, dt.day)
                log.fatal(
                    f"Frequency {freq} antenna corrections for satellite PRN {ant} valid from {valid_from} is not unique in file {self.file_path}."
                )

            self.data[ant][dt].update(tmp)

        # Save receiver antenna correction in data structure
        else:
            if freq in self.data[ant]:
                log.fatal(
                    f"Frequency {freq} antenna corrections for receiver antenna {ant} is not unique in file {self.file_path}."
                )

            self.data[ant].update(tmp)

        cache["num_freq_counter"] = cache["num_freq_counter"] + 1

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return []
