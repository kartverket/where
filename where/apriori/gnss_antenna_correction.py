"""Handling of GNSS antenna corrections

Description:
------------

The module includes a class for handling apriori GNSS antenna corrections. Read GNSS ANTEX file (see
:cite:`hilla2010`).

"""
# Standard library imports
import datetime
from typing import Dict, Optional, Union

# External library imports
from collections import UserDict
import numpy as np

# Migard imports
from midgard.collections import enums
from midgard.dev import plugins

# Where imports
from where import parsers
from where.lib import log
from where.data import position


@plugins.register
class AntennaCorrection(UserDict):
    """A class for representing apriori GNSS antenna correction data

    The attribute "data" is a dictionary with GNSS satellite PRN or receiver antenna as key. The GNSS satellite antenna
    corrections are time dependent and saved with "valid from" datetime object entry. The dictionary looks like:

        dout = { <prn> : { <valid from>: { cospar_id:   <value>,
                                           sat_code:    <value>,
                                           sat_type:    <value>,
                                           valid_until: <value>,
                                           azimuth:     <list with azimuth values>,
                                           elevation:   <list with elevation values>,
                                           <frequency>: { azi: [<list with azimuth-elevation dependent corrections>],
                                                          neu: [north, east, up],
                                                          noazi: [<list with elevation dependent corrections>] }}},

                 <receiver antenna> : { azimuth:     <list with azimuth values>,
                                        elevation:   <list with elevation values>,
                                        <frequency>: { azi: [<array with azimuth-elevation dependent corrections>],
                                                       neu: [north, east, up],
                                                       noazi: [<list with elevation dependent corrections>] }}}

    with following entries:

    | Value              | Type              | Description                                                             |
    | :----------------- | :---------------- | :---------------------------------------------------------------------- |
    | azi                | numpy.ndarray     | Array with azimuth-elevation dependent antenna correction in [mm] with  |
    |                    |                   | the shape: number of azimuth values x number of elevation values.       |
    | azimuth            | numpy.ndarray     | List with azimuth values in [rad] corresponding to antenna corrections  |
    |                    |                   | given in `azi`.                                                         |
    | cospar_id          | str               | COSPAR ID <yyyy-xxxa>: yyyy -> year when the satellite was put in       |
    |                    |                   | orbit, xxx -> sequential satellite number for that year, a -> alpha     |
    |                    |                   | numeric sequence number within a launch                                 |
    | elevation          | numpy.ndarray     | List with elevation values in [rad] corresponding to antenna            |
    |                    |                   | corrections given in `azi` or `noazi`.                                  |
    | <frequency>        | str               | Frequency identifier (e.g. G01 - GPS L1)                                |
    | neu                | list              | North, East and Up eccentricities in [m]. The eccentricities of the     |
    |                    |                   | mean antenna phase center is given relative to the antenna reference    |
    |                    |                   | point (ARP) for receiver antennas or to the center of mass of the       |
    |                    |                   | satellite in X-, Y- and Z-direction.                                    |
    | noazi              | numpy.ndarray     | List with elevation dependent (non-azimuth-dependent) antenna           |
    |                    |                   | correction in [mm].                                                     |
    | <prn>              | str               | Satellite code e.g. GPS PRN, GLONASS slot or Galileo SVID number        |
    | <receiver antenna> | str               | Receiver antenna name together with radome code                         |
    | sat_code           | str               | Satellite code e.g. GPS SVN, GLONASS number or Galileo GSAT number      |
    | sat_type           | str               | Satellite type (e.g. BLOCK IIA)                                         |
    | valid_from         | datetime.datetime | Start of validity period of satellite in GPS time                       |
    | valid_until        | datetime.datetime | End of validity period of satellite in GPS time                         |


    Attributes:
        data (dict):           Data read from GNSS Antenna Exchange (ANTEX) file
        file_path (str):       ANTEX file path

    Methods:
        satellite_phase_center_offset(): Determine satellite phase center offset correction vectors given in ITRS
        satellite_type(): Get satellite type from ANTEX file (e.g. BLOCK IIF, GALILEO-1, GALILEO-2, GLONASS-M,
                          BEIDOU-2G, ...)
        _used_date(): Choose correct date for use of satellite antenna corrections
    """

    def __init__(self, file_key: Optional[str] = "antex") -> None:
        """Set up a new GNSS antenna correction object by parsing ANTEX file

        The parsing is done by :mod:`where.parsers.antex`.
        """
        parser = parsers.parse_key_existing(file_key)
        self.data = parser.as_dict()
        self.file_path = parser.file_path
        log.info(f"Read ANTEX file: {self.file_path}")


    def satellite_phase_center_offset(
        self, dset: "Dataset", sys_freq: Union[None, Dict[str, Dict[str, str]]] = None
    ) -> "PosVelDeltaArray":
        """Determine satellite phase center offset correction vectors given in ITRS

        Satellite phase center offset (PCO) corrections are frequency dependent. The argument 'sys_freq' defines, which
        frequencies for a given GNSS should be used. If two frequencies for a GNSS are given, then the PCO is
        determined as ionospheric-free linear combination. If 'sys_freq' is not defined as input argument, then 
        'sys_freq' is generated based on the given observation types in dataset 'dset'.

        Args:
            dset:      Model data.
            sys_freq:  Dictionary with frequency or frequency combination given for GNSS identifier:
                         sys_freq = { <sys_id>: <freq> }  (e.g. sys_freq = {'E': 'E1',  'G': 'L1_L2'} )

        Returns:
            Satellite phase center offset correction vectors given in ITRS in meter

        """
        # GNSS          Freq number      GNSS freq
        #               L<num>/C<num>
        # ___________________________________________
        # C (BeiDou):   2                'B1'
        #               7                'B2'
        #               6                'B3'
        # G (GPS):      1                'L1'
        #               2                'L2'
        #               5                'L5'
        # R (GLONASS):  1                'G1'
        #               2                'G2'
        #               3                'G3'
        # E (Galileo):  1                'E1'
        #               8                'E5 (E5a+E5b)'
        #               5                'E5a'
        #               7                'E5b'
        #               6                'E6'
        # I (IRNSS):    5                'L5'
        #               9                'S'
        # J (QZSS):     1                'L1'
        #               2                'L2'
        #               5                'L5'
        #               6                'LEX'
        # S (SBAS):     1                'L1'
        #               5                'L5'
        obstype_to_gnss_freq = {
            "C": {"2": "B1", "7": "B2", "6": "B3"},
            "E": {"1": "E1", "8": "E5", "5": "E5a", "7": "E5b", "6": "E6"},
            "G": {"1": "L1", "2": "L2", "5": "L5"},
            "I": {"5": "L5", "9": "S"},
            "J": {"1": "L1", "2": "L2", "5": "L5", "6": "LEX"},
            "R": {"1": "G1", "2": "G2", "3": "G3"},
            "S": {"1": "L1", "5": "L5"},
        }

        correction = position.PosVelDelta(
            val=np.zeros((dset.num_obs, 6)), system="yaw", ref_pos=dset.sat_posvel, time=dset.time
        )
        used_date = None

        # Get GNSS frequency based on observation type
        if sys_freq is None:
            sys_freq = dict()
            obstypes = dset.meta["obstypes"]
            for sys in obstypes:
                freqs = {o[1] for o in obstypes[sys]}
                sys_freq[sys] = "_".join(f"{obstype_to_gnss_freq[sys][f]}" for f in sorted(freqs))

        # Loop over all satellites given in RINEX observation file and configuration file
        for sat in dset.unique("satellite"):

            # Skip satellites, which are not given in ANTEX file
            if sat not in self.data:
                log.warn(
                    f"Satellite {sat} is not given in ANTEX file {self.file_path}. That means no satellite antenna "
                    f"phase center offset correction can be applied for satellite {sat}."
                )
                continue

            # Get array with information about, when observation are available for the given satellite (indicated by True)
            idx = dset.filter(satellite=sat)

            # Get used date
            used_date = self._used_date(sat, dset.analysis["rundate"])
            if used_date is None:
                continue

            # Add PCO to Dataset meta
            pco_sat = self.get_pco_sat(sat, sys_freq, used_date)
            dset.meta.setdefault("pco_sat", dict()).update({sat: pco_sat.tolist()})
            correction[idx] = np.repeat(np.append(pco_sat, np.zeros(3))[None, :], np.sum(idx), axis=0)

        # Transform PCO given in satellite body-fixed reference frame (for GPS and Galileo assumed to be aligned
        # to yaw-steering reference frame) to ITRS
        return correction.trs


    def get_pco_sat(self, sat: str, sys_freq: Dict[str, Dict[str, str]], used_date: datetime.datetime) -> np.ndarray:
        """Get satellite PCO in satellite reference system

        If two frequencies are given over the 'sys_freq' argument, then the PCOs are determined as an ionospheric linear
        combination.

        Args:
            sat:        Satellite identifier.
            sys_freq:   Dictionary with frequency or frequency combination given for GNSS identifier:
                            sys_freq = { <sys_id>: <freq> }  (e.g. sys_freq = {'E': 'E1',  'G': 'L1_L2'} )

            used_date:  Correct date for use of satellite antenna corrections

        Returns:
            Satellite PCO in satellite reference system
        """
        # GNSS          GNSS freq        ANTEX freq
        # ___________________________________________
        # C (BeiDou):   B1               'C02'
        #               B2               'C07'
        #               B3               'C06'
        # G (GPS):      L1               'G01'
        #               L2               'G02'
        #               L5               'G05'
        # R (GLONASS):  G1               'R01'
        #               G2               'R02'
        # E (Galileo):  E1               'E01'
        #               E5 (E5a+E5b)     'E08'
        #               E5a              'E05'
        #               E5b              'E07'
        #               E6               'E06'
        # I (IRNSS):    L5               'I05'
        #               S                'I09'
        # J (QZSS):     L1               'J01'
        #               L2               'J02'
        #               L5               'J05'
        #               LEX              'J06'
        # S (SBAS):     L1               'S01'
        #               L5               'S05'
        gnss_to_antex_freq = {
            "C": {"B1": "C02", "B2": "C07", "B3": "C06"},
            "E": {"E1": "E01", "E5": "E08", "E5a": "E05", "E5b": "E07", "E6": "E06"},
            "G": {"L1": "G01", "L2": "G02", "L5": "G05"},
            "I": {"L5": "I05", "S": "I09"},
            "J": {"L1": "J01", "L2": "J02", "L5": "J05", "LEX": "J06"},
            "R": {"G1": "R01", "G2": "R02", "G3": "R03"},
            "S": {"L1": "S01", "L5": "S05"},
        }

        # Get GNSS and frequency/frequencies
        sys = sat[0]  # GNSS identifier
        freq = sys_freq[sys]
        if "_" in freq:
            freq = freq.split("_")
        else:
            freq = [freq]

        # Get satellite PCO for one frequency
        if len(freq) == 1:

            # Get satellite phase center offset (PCO) given in satellite reference system
            pco_sat = np.array(self.data[sat][used_date][gnss_to_antex_freq[sys][freq[0]]]["neu"])

            log.debug(f"PCO of satellite {sat} for frequency {sys}:{sys_freq[sys]}: {pco_sat.tolist()[0]}.")

        # Get satellite PCO for ionospheric-free linear combination based on two-frequencies
        elif len(freq) == 2:

            # Coefficient of ionospheric-free linear combination
            f1 = getattr(enums, "gnss_freq_" + sys)[freq[0]]  # Frequency of 1st band
            f2 = getattr(enums, "gnss_freq_" + sys)[freq[1]]  # Frequency of 2nd band
            n = f1 ** 2 / (f1 ** 2 - f2 ** 2)
            m = -f2 ** 2 / (f1 ** 2 - f2 ** 2)

            # Get satellite phase center offset (PCO) given in satellite reference system
            pco_sat_f1 = np.array(self.data[sat][used_date][gnss_to_antex_freq[sys][freq[0]]]["neu"])
            pco_sat_f2 = np.array(self.data[sat][used_date][gnss_to_antex_freq[sys][freq[1]]]["neu"])

            # Generate ionospheric-free linear combination
            pco_sat = n * pco_sat_f1 + m * pco_sat_f2

            log.debug(
                f"Ionospheric-free linear combination PCOs of satellite {sat} for frequency combination "
                f"{sys}:{sys_freq[sys]}:  {pco_sat.tolist()[0]}."
            )

        else:
            log.fatal(
                f"Wrong frequency type '{sys}:{'_'.join(freq)}'. Note: 'signals' configuration option should represent "
                f"one- or two-frequencies (e.g. E:E1 or E:E1_E5a)."
            )

        return pco_sat


    def get_satellite_info(self, satellite: str, date: datetime.datetime) -> Dict[str, str]:
        """Get satellite information for a given date

        Args:
            sat:   Satellite identifier.
            date:  Date.

        Returns:
            Satellite information
        """
        if satellite not in self.data.keys():
            log.fatal(
                f"Satellite '{satellite}' is not given in ANTEX file {self.file_path}."
            )

        used_date = self._used_date(satellite, date)

        return {
            "sat_type": self.data[satellite][used_date]["sat_type"],
            "sat_code": self.data[satellite][used_date]["sat_code"],
            "cospar_id": self.data[satellite][used_date]["cospar_id"],
        }


    def satellite_type(self, dset: "Dataset") -> np.core.defchararray.chararray:
        """Get satellite type from ANTEX file (e.g. BLOCK IIF, GALILEO-1, GALILEO-2, GLONASS-M, BEIDOU-2G, ...)

        Args:
            dset:   Model data.

        Returns:
            Satellite type information

        """
        sat_types = np.zeros(dset.num_obs, dtype=object)
        used_date = None

        # Loop over all satellites given in RINEX observation file and configuration file
        for sat in dset.unique("satellite"):

            # Skip satellites, which are not given in ANTEX file
            if sat not in self.data:
                log.warn(
                    f"Satellite {sat} is not given in ANTEX file {self.file_path}. That means no satellite antenna "
                    f"phase center offset correction can be applied for satellite {sat}."
                )
                continue

            # Get array with information about, when observation are available for the given satellite (indicated
            # by True)
            idx = dset.filter(satellite=sat)

            # Get used date
            used_date = self._used_date(sat, dset.analysis["rundate"])
            if used_date is None:
                continue

            # Get satellite type
            sat_type = self.data[sat][used_date]["sat_type"]
            sat_types[idx] = sat_type

        return sat_types.astype(np.str_)


    def _used_date(self, satellite: str, given_date: datetime.date) -> Union[datetime.datetime, None]:
        """Choose correct date for use of satellite antenna corrections

        Satellite antenna correction are time dependent.

        Args:
            satellite (str):              Satellite identifier.
            given_date (datetime.date):   Given date used for finding corresponding time period in ANTEX file

        Returns:
            Date for getting correct satellite antenna corrections related to given date
        """
        used_date = None
        # TODO: Would it be not better to define rundate as datetime.datetime instead datetime.date?
        given_date = datetime.datetime.combine(given_date, datetime.time())  # conversion from date to datetime

        for date in sorted(self.data[satellite]):
            if date <= given_date:
                used_date = date

        if (used_date is None) or (given_date > self.data[satellite][used_date]["valid_until"]):
            log.warn(f"No satellite phase center offset is given for satellite {satellite} and date {given_date}.")

        return used_date
