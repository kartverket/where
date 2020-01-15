"""Get apriori data for GNSS satellite clocks

Description:
------------

GNSS satellite clock data are read from RINEX clock file.




"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import data
from where import parsers


@plugins.register
def get_gnss_sat_clk(rundate):
    """Get GNSS satellite clock data

    The parsing is done by :mod:`where.parsers.rinex_clk_parser`.

    Returns:
          Dataset: Dataset with precise satellite clock values stored as follows

    ====================  ===============  =======  ========================================================
     Field                 Type             Unit     Description
    ====================  ===============  =======  ========================================================
     sat_clock_bias        numpy.ndarray     m       Satellite clock offset from GPS time
     satellite             list                      Satellite PRN number
     time                  TimeTable                 Observation epoch in GPS time
     time.data             Time                      TODO
     time.utc              Time                      TODO
    ====================  ===============  =======  ========================================================

    """
    gnss_sat_clk_data = data.Dataset(
        rundate, tech=None, stage=None, dataset_name="gnss_sat_clk", dataset_id=0, empty=True
    )
    parser = parsers.parse("rinex_clk", rundate=rundate)
    parser.write_to_dataset(gnss_sat_clk_data)

    return gnss_sat_clk_data
