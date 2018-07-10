"""Handling of GNSS antenna corrections

Description:
------------
The module includes a class for handling apriori GNSS antenna corrections. Read GNSS ANTEX file (see :cite:`hilla2010`).


"""
# Standard library imports
import datetime

# External library imports
from collections import UserDict
import numpy as np

# Where imports
from where import data
from where import parsers
from where.lib import constant
from where.lib import log
from where.lib import plugins


@plugins.register
class AntennaCorrection(UserDict):
    """A class for representing apriori GNSS antenna correction data

    The attribute "data" is a dictionary with GNSS satellite PRN or receiver antenna as key. The GNSS satellite antenna
    corrections are time dependent and saved with "valid from" datetime object entry. The dictionary looks like:

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


    Attributes:
        data (dict):            Data read from GNSS Antenna Exchange (ANTEX) file
        file_path (str):       ANTEX file path

    Methods:
        satellite_phase_center_offset(): Determine satellite phase center offset correction vectors given in ITRS
        satellite_type(): Get satellite type from ANTEX file (e.g. BLOCK IIF, GALILEO-1, GALILEO-2, GLONASS-M, 
                          BEIDOU-2G, ...)
        _used_date(): Choose correct date for use of satellite antenna corrections
    """

    def __init__(self, file_key="gnss_antex"):
        """Set up a new GNSS antenna correction object by parsing ANTEX file

        The parsing is done by :mod:`where.parsers.gnss_antex`.
        """
        parser = parsers.parse_key(file_key)
        self.data = parser.as_dict()
        self.file_path = parser.file_path

    def satellite_phase_center_offset(self, dset, sys_freq):
        """Determine satellite phase center offset correction vectors given in ITRS

        Satellite phase center offset (PCO) corrections are frequency dependent. The argument 'sys_freq' defines, which
        frequencies for a given GNSS should be used. If two frequencies for a GNSS are given, then the PCO is
        determined as ionospheric-free linear combination.

        Args:
            dset (Dataset):   Model data.
            sys_freq (dict):  Dictionary with frequency or frequency combination given for GNSS
                              identifier:
                                sys_freq = { <sys_id>: <freq> }  (e.g. sys_freq = {'E': 'E1',  'G': 'L1_L2'} )

        Returns:
            numpy.ndarray: Satellite phase center offset correction vectors given in ITRS in meter

        """
        correction = np.zeros((dset.num_obs, 3))
        used_date = None

        # Loop over all satellites given in RINEX observation file and configuration file
        for sat in dset.unique("satellite"):

            # Skip satellites, which are not given in ANTEX file
            if sat not in self.data:  # antex:
                log.warn(
                    "Satellite {} is not given in ANTEX file {}. That means no satellite antenna phase center offset "
                    "correction can be applied for satellite {}.",
                    sat,
                    self.file_path,
                    sat,
                )
                continue

            # Get array with information about, when observation are available for the given satellite (indicated by True)
            idx = dset.filter(satellite=sat)

            # Get used date
            used_date = self._used_date(sat, dset.rundate)
            if used_date is None:
                continue

            # Add PCO to Dataset meta
            system = sat[0]
            pco_sat = self._get_pco_sat(sat, sys_freq, used_date)
            dset.meta.setdefault("pco_sat", dict()).update({sat: pco_sat.tolist()[0]})

            # Transform PCO given in satellite body-fixed reference frame (for GPS and Galileo assumed to be aligned
            # to yaw-steering reference frame) to ITRS
            pco_itrs = dset.sat_posvel.convert_yaw_to_itrs(pco_sat)
            # pco_itrs = dset.sat_posvel._yaw2itrs[idx][0] @ np.array(pco_sat)

            correction[idx] = pco_itrs[idx]

        return correction

    def _get_pco_sat(self, sat, sys_freq, used_date):
        """Get satellite PCO in satellite reference system

        If two frequencies are given over the 'sys_freq' argument, then the PCOs are determined as an ionospheric linear
        combination.

        Args:
            sat (str):                      Satellite identifier.
            sys_freq (dict):                Dictionary with frequency or frequency combination given for GNSS
                                            identifier:
                                                sys_freq = { <sys_id>: <freq> }  
                                                (e.g. sys_freq = {'E': 'E1',  'G': 'L1_L2'} )
            used_date (datetime.datetime):  Correct date for use of satellite antenna corrections     

        Returns:
            numpy.ndarray: Satellite PCO in satellite reference system
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
            "R": {"G1": "R01", "G2": "R02"},
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
            pco_sat = np.array([self.data[sat][used_date][gnss_to_antex_freq[sys][freq[0]]]["neu"]])

        # Get satellite PCO for ionospheric-free linear combination based on two-frequencies
        elif len(freq) == 2:

            # Coefficient of ionospheric-free linear combination
            f1 = constant.get("gnss_freq_" + sys, source=freq[0])  # Frequency of 1st band
            f2 = constant.get("gnss_freq_" + sys, source=freq[1])  # Frequency of 2nd band
            n = f1 ** 2 / (f1 ** 2 - f2 ** 2)
            m = -f2 ** 2 / (f1 ** 2 - f2 ** 2)

            # Get satellite phase center offset (PCO) given in satellite reference system
            pco_sat_f1 = np.array([self.data[sat][used_date][gnss_to_antex_freq[sys][freq[0]]]["neu"]])
            pco_sat_f2 = np.array([self.data[sat][used_date][gnss_to_antex_freq[sys][freq[1]]]["neu"]])

            # Generate ionospheric-free linear combination
            pco_sat = n * pco_sat_f1 + m * pco_sat_f2

        else:
            log.fatal(
                f"Wrong frequency type '{sys}:{'_'.join(freq)}'. Note: 'signals' configuration option should represent one- or two-frequencies (e.g. E:E1 or E:E1_E5a)."
            )

        return pco_sat

    def satellite_type(self, dset):
        """Get satellite type from ANTEX file (e.g. BLOCK IIF, GALILEO-1, GALILEO-2, GLONASS-M, BEIDOU-2G, ...)

        Args:
            dset (Dataset):   Model data.

        Returns:
            numpy.core.defchararray.chararray: Satellite type information

        """
        sat_types = data.text_table.TextList([""]) * dset.num_obs
        used_date = None

        # Loop over all satellites given in RINEX observation file and configuration file
        for sat in dset.unique("satellite"):

            # Skip satellites, which are not given in ANTEX file
            if sat not in self.data:
                log.warn(
                    "Satellite {} is not given in ANTEX file {}. That means no satellite antenna phase center offset "
                    "correction can be applied for satellite {}.",
                    sat,
                    self.file_path,
                    sat,
                )
                continue

            # Get array with information about, when observation are available for the given satellite (indicated by True)
            idx = dset.filter(satellite=sat)

            # Get used date
            used_date = self._used_date(sat, dset.rundate)
            if used_date is None:
                continue

            # Get satellite type
            sat_type = self.data[sat][used_date]["sat_type"]
            sat_types[idx] = sat_type

        return sat_types

    def _used_date(self, satellite, given_date):
        """Choose correct date for use of satellite antenna corrections

        Satellite antenna correction are time dependent.

        Args:
            satellite (str):              Satellite identifier.
            given_date (datetime.date):   Given date used for finding corresponding time period in ANTEX file
        """
        used_date = None
        # TODO: Would it be not better to define rundate as datetime.datetime instead datetime.date?
        given_date = datetime.datetime.combine(given_date, datetime.time())  # conversion from date to datetime

        for date in sorted(self.data[satellite]):
            if date <= given_date:
                used_date = date

        if (used_date is None) or (given_date > self.data[satellite][used_date]["valid_until"]):
            log.warn(
                "No satellite phase center offset is given for satellite {} and date {}.",
                satellite,
                given_date,
                satellite,
            )

        return used_date
