"""A parser for reading GNSS raw data from `GnssLogger` Android App

Example:
--------

    from where import parsers
    parser = parsers.parse('gnss_android_raw_data', rundate=rundate)

Description:
------------

Reads raw data file from `GnssLogger` Android App.

"""

# Standard library imports
import itertools

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers import parser
from midgard.math.constant import constant
from where.lib import gnss
from where.lib.unit import Unit


@plugins.register
class GnssAndroidRawDataParser(parser.Parser):
    """
    """

    def __init__(self, rundate, station=None):
        """
        Args:
            rundate:       The model run date.
        """
        super().__init__(rundate)
        self.file_key = "gnss_android_raw_data"

        if station:
            self.vars["station"] = station.lower()
            self.vars["STATION"] = station.upper()

    #
    # PARSERS
    #
    def setup_parsers(self):
        """Parsers defined for reading GNSS raw data file line by line.
        """
        file_parser = parser.define_parser(
            end_marker=lambda _l, _ln, _n: True,
            label=lambda line, _ln: line[0:3].strip(),
            parser_def={
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8
                # Fix,gps,37.422541,-122.081659,-33.000000,0.000000,3.000000,1467321969000
                "Fix": {
                    "parser": self.parse_fix,
                    "fields": [
                        "dummy",
                        "provider",
                        "latitude",
                        "longitude",
                        "altitude",
                        "speed",
                        "accuracy",
                        "time_utc",
                    ],
                },
                # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0
                # Raw,72065126,72076939000000,,,-1151285108458178048,0.0,26.542398700257763,-0.634638974724185, ...
                "Raw": {
                    "parser": self.parse_raw,
                    "delimiter": ",",
                    "fields": [
                        "dummy",
                        "ElapsedRealtimeMillis",
                        "TimeNanos",
                        "LeapSecond",
                        "TimeUncertaintyNanos",
                        "FullBiasNanos",
                        "BiasNanos",
                        "BiasUncertaintyNanos",
                        "DriftNanosPerSecond",
                        "DriftUncertaintyNanosPerSecond",
                        "HardwareClockDiscontinuityCount",
                        "Svid",
                        "TimeOffsetNanos",
                        "State",
                        "ReceivedSvTimeNanos",
                        "ReceivedSvTimeUncertaintyNanos",
                        "Cn0DbHz",
                        "PseudorangeRateMetersPerSecond",
                        "PseudorangeRateUncertaintyMetersPerSecond",
                        "AccumulatedDeltaRangeState",
                        "AccumulatedDeltaRangeMeters",
                        "AccumulatedDeltaRangeUncertaintyMeters",
                        "CarrierFrequencyHz",
                        "CarrierCycles",
                        "CarrierPhase",
                        "CarrierPhaseUncertainty",
                        "MultipathIndicator",
                        "SnrInDb",
                        "ConstellationType",
                    ],
                },
            },
        )

        return itertools.repeat(file_parser)

    def parse_fix(self, line, _):
        """Parse 'Fix' entries of GNSS raw data file to instance variable 'data'.
        """
        for k, v in line.items():
            if k is "dummy":
                continue
            if k is "provider":
                self.data.setdefault(k, list()).append(v)
                continue

            self.data.setdefault(k, list()).append(float(v))

    def parse_raw(self, line, _):
        """Parse 'Raw' entries of GNSS raw data file to instance variable 'data'.
        """
        # TODO: Make parsing depending of 'State' value. What do the 'State' numbers mean?
        constellationType = {"0": None, "1": "G", "2": "S", "3": "R", "4": "J", "5": "C", "6": "E"}

        for k, v in line.items():
            if k is "dummy":
                continue
            if k is "Svid":
                system = constellationType[line["ConstellationType"]]
                if system == None:
                    log.warn("GNSS is unknown.")
                    continue
                self.data.setdefault("system", list()).append(system)
                self.data.setdefault("satellite", list()).append(system + str(v).zfill(2))
            if v is "":
                v = float("nan")

            self.data.setdefault(k, list()).append(float(v))

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return [self.determine_pseudorange, self.get_sitepos]

    def determine_pseudorange(self):
        """Determine pseudorange based on ION 2016 tutorial "Raw GNSS Measurements from Android Phones".
        """

        # Determine GPS week
        week = np.floor(-np.array(self.data["FullBiasNanos"]) * Unit.nanosecond2second / 604800)

        # GNSS signal arriving time at measurement time (GPS time) referenced to GPS week
        tRxNanos = (
            (np.array(self.data["TimeNanos"], dtype=float) + np.array(self.data["TimeOffsetNanos"], dtype=float))
            - (np.array(self.data["FullBiasNanos"], dtype=float) + np.array(self.data["BiasNanos"], dtype=float))
            - (week * 604800e9)
        )

        if np.all(tRxNanos >= 604800e9):
            log.fatal("tRxNanos should be <= GPS nanoseconds.")
        if np.all(tRxNanos <= 0.0):
            log.fatal("tRxNanos should be >= 0.")

        self.data["week"] = week
        self.data["tRxNanos"] = tRxNanos
        self.data["time"] = gnss.gpssec2jd(week, tRxNanos * Unit.ns2s)  # TODO: Better solution -> time.py?

        # GNSS satellite transmission time at measurement time (GPS time) referenced to GPS week
        tTxNanos = np.array(self.data["ReceivedSvTimeNanos"], dtype=float)
        self.data["sat_time"] = gnss.gpssec2jd(week, tTxNanos * Unit.nanosecond2second)
        # TODO: Check GPS week rollover (see ProcessGnssMeas.m)

        self.data["pseudorange"] = (tRxNanos - tTxNanos) * Unit.nanosecond2second * constant.c  # in meters

    def get_sitepos(self):
        """Determine site position by converting given longitude, latitude and height for a station to geocentric
           cartesian coordinates
        """
        if "latitude" in self.data.keys():
            # TODO: Better solution?
            x, y, z = gnss.llh2xyz(
                np.deg2rad(self.data["latitude"][0]), np.deg2rad(self.data["longitude"][0]), self.data["altitude"][0]
            )
        else:
            x = y = z = 0.0

        for idx in range(0, len(self.data["time"])):
            self.data.setdefault("site_pos", list()).append([x, y, z])

    #
    # WRITE DATA
    #
    def write_to_dataset(self, dset):
        """Store GNSS data in a dataset

        Args:
            dset: The Dataset where data are stored.
        """
        dset.num_obs = len(self.data["pseudorange"])

        dset.add_float("C1C", val=self.data["pseudorange"])

        dset.add_time("time", val=self.data["time"], scale="gps", format="jd")
        # dset.add_time('sat_time', val=self.data['sat_time'], scale='gps', format='jd')

        dset.add_text("station", val=np.full(dset.num_obs, "android", dtype=str))
        dset.add_text("satellite", val=self.data["satellite"])
        dset.add_text("satnum", val=self.data["Svid"])
        dset.add_text("system", val=self.data["system"])

        dset.add_position("site_pos", time="time", itrs=np.array(self.data["site_pos"]))
