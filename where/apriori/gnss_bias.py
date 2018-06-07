"""Handling of GNSS biases

Description:
------------
The module includes a class for handling apriori GNSS bias corrections. Read SINEX bias file (see :cite:`sinex_bias`).


"""
# Standard library imports
import datetime

# External library imports
from collections import UserDict
import numpy as np

# Where imports
from where import data
from where import parsers
from where.lib import config
from where.lib import log
from where.lib import plugins


@plugins.register
class BiasCorrection(UserDict):
    """A class for representing apriori GNSS bias correction data

    The attribute 'data' is a dictionary with satellite PRNs as key and a dictionary with GNSS code bias combinations
    as value. The GNSS code bias combinations are time dependent and saved with 'start_time' datetime object entry. The
    dictionary looks like:

        data = { '<prn>': {'<code1>-<code2>': {<start_time>: {end_time: <date>,
                                                              estimate: <dsb_estimate>,
                                                              sigma: <dsb_sigma>}}}
                           svn: <number>}

    Example:
                                              
        data = { 'G29': {'C1C-C1W': {datetime.datetime(2018, 2, 1, 0, 0): {'end_time': '2018:033:00000',
                                                                          'estimate': -5.410000000000001e-10,
                                                                          'sigma': 8.500000000000001e-12}},
                        'C1C-C2W': {datetime.datetime(2018, 2, 1, 0, 0): {'end_time': '2018:033:00000',
                                                                          'estimate': 1.4690000000000002e-09,
                                                                          'sigma': 2.0500000000000004e-11}},
                        ....
                        'svn': 'G057'},
                'G30': {'C1C-C1W': {datetime.datetime(2018, 2, 1, 0, 0): {'end_time': '2018:033:00000',
                                                                          'estimate': 4.69e-10,
                                                                          'sigma': 8.500000000000001e-12}},
                ...}

    with following entries:

    =================== =================== ========================================================================
     Value               Type                Description
    =================== =================== ========================================================================
     <code1>-<code2>     str                 GNSS code bias combination (e.g. 'C1C-C1W')
     estimate            float               Differential Signal Bias (DSB) estimate in [s]
     <prn>               str                 Satellite code e.g. GPS PRN, GLONASS slot or Galileo SVID number
     <start_time>        datetime.datetime   Start of validity period of GNSS code bias
     sigma               float               Differential Signal Bias (DSB) estimated standard deviation in [s]
     svn                 str                 Satellite SVN code.
     end_time            datetime.datetime   End of validity period of GNSS code bias
    =================== =================== =======================================================================


    Attributes:
        data (dict):            Data read from GNSS bias (SINEX) file
        file_path (str):        GNSS SINEX bias file path
    """

    def __init__(self, rundate, file_key="gnss_sinex_bias"):
        """Set up a new GNSS bias object by parsing SINEX bias file

        The parsing is done by :mod:`where.parsers.gnss_sinex_bias`.
        """
        parser = parsers.parse_key(file_key=file_key, file_vars=config.date_vars(rundate))
        self.data = parser.as_dict()

    def get_dsb(self, satellite, dsb, given_date):
        """Get correct Differential Signal Bias for given satellite, DSB code and date

        Args:
            satellite (str):              Satellite identifier.
            dsb (str):                    Differential Signal Bias combination (e.g. 'C1C-C1W')
            given_date (datetime.date):   Given date used for finding corresponding time period in GNSS bias file
        """
        used_date = None
        # TODO: Would it be not better to define rundate as datetime.datetime instead datetime.date?
        given_date = datetime.datetime.combine(given_date, datetime.time())  # conversion from date to datetime

        for date in sorted(self.data[satellite][dsb]):
            if date <= given_date:
                used_date = date

        if (used_date is None) or (given_date > self.data[satellite][dsb][used_date]["end_time"]):
            log.warn(
                "No satellite GNSS bias is given for satellite {}, Differential Signal Bias (DSB) {} and time {}. ",
                satellite,
                dsb,
                given_date,
            )

        return self.data[satellite][dsb][used_date]
