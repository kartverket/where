"""A parser for reading SINEX bias file for GNSS biases

Example:
--------

    from where import parsers
    p = parsers.parse_file(parser_name='gnss_sinex_bias', file_path='CAS0MGXRAP_20180320000_01D_01D_DCB.BSX.gz')
    data = p.as_dict()

Description:
------------

Reads data from files in SINEX GNSS bias file format 1.0 (see :cite:`sinex_bias`).


"""
# Standard library imports
from datetime import datetime, timedelta
import itertools

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.parsers._parser_chain import ChainParser, ParserDef
from where.lib import time
from where.lib.unit import Unit


# TODO: Should we use SINEX class instead?
@plugins.register
class GnssBiasParser(ChainParser):
    """A parser for reading SINEX bias file for GNSS biases

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
    =================== =================== ========================================================================
    """

    #
    # PARSER
    #
    def setup_parser(self):
        """Parsers defined for reading SINEX bias file line by line
        """
        #
        # +BIAS/SOLUTION
        # ----+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1
        # *BIAS SVN_ PRN STATION__ OBS1 OBS2 BIAS_START____ BIAS_END______ UNIT __ESTIMATED_VALUE____ _STD_DEV___
        #  DSB  G063 G01           C1C  C1W  2018:032:00000 2018:033:00000 ns                 -0.9800      0.0085
        #  DSB  G061 G02           C1C  C1W  2018:032:00000 2018:033:00000 ns                  1.3100      0.0085
        #  DSB  G069 G03           C1C  C1W  2018:032:00000 2018:033:00000 ns                 -1.4730      0.0090

        # TODO: How to handle different SINEX BIAS format versions?
        file_parser = ParserDef(
            # TODO: start_marker=lambda line, _ln: startswith('+BIAS/SOLUTION') ????
            end_marker=lambda line, _ln, _n: line.startswith("-BIAS/SOLUTION"),
            label=lambda line, _ln: line.startswith(" DSB"),
            parser_def={
                True: {
                    "parser": self.parse_bias_solution,
                    "fields": {
                        "bias": (0, 5),
                        "svn": (5, 10),
                        "prn": (10, 14),
                        "station": (14, 25),  # TODO: Old format is not consistent with SINEX BIAS format 1.0.
                        "obs1": (25, 29),
                        "obs2": (29, 34),
                        "start_time": (
                            34,
                            49,
                        ),  # TODO: Old format is not consistent with SINEX BIAS format 1.0. -> yy instead of yyyy
                        "end_time": (
                            49,
                            64,
                        ),  # TODO: Old format is not consistent with SINEX BIAS format 1.0. -> yy instead of yyyy
                        "unit": (64, 69),
                        "estimate": (69, 91),
                        "estimate_sigma": (91, 103),
                        "slope_estimate": (103, 125),  # TODO: Not implemented.
                        "slope_estimate_sigma": (125, 137),  # TODO: Not implemented.
                    },
                }
            },
        )

        return itertools.repeat(file_parser)

    def parse_bias_solution(self, line, _):
        """Parse BIAS/SOLUTION block entries to instance variable 'data'
        """
        # Prepare data dictionary entries
        start_time = datetime.strptime(line["start_time"][0:8], "%Y:%j") + timedelta(
            seconds=int(line["start_time"][9:14])
        )
        end_time = datetime.strptime(line["end_time"][0:8], "%Y:%j") + timedelta(seconds=int(line["end_time"][9:14]))

        if line["unit"] == "cyc":  # TODO: Differential phase bias conversion is not handled yet in Where.
            print("WARN: Differential phase bias conversion is not handled yet in Where.")
            return 0  # TODO: Is it correct to handle it like that.
        elif line["unit"] == "ns":
            estimate = float(line["estimate"]) * Unit.nanosecond2second
            sigma = float(line["estimate_sigma"]) * Unit.nanosecond2second
        else:
            print(
                "FATAL: Unit '{}' in file {} is not defined in SINEX bias format."
                "".format(line["unit"], self.file_path)
            )
            return 0  # TODO: How should it be handled? sys.exit(1)?

        # Update data dictionary with current DSB line
        self.data.setdefault(line["prn"], dict())
        self.data[line["prn"]].update({"svn": line["svn"]})
        self.data[line["prn"]].setdefault(line["obs1"] + "-" + line["obs2"], dict())
        self.data[line["prn"]][line["obs1"] + "-" + line["obs2"]].setdefault(start_time, dict()).update(
            {"estimate": estimate, "sigma": sigma, "end_time": end_time}
        )

    #
    # SETUP CALCULATION
    #
    def setup_calculators(self):
        """List steps necessary for postprocessing
        """
        return []
