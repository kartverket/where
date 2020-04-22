"""Get apriori data for gravitational deformation of VLBI antennas

Description:
    Reads the gravitational deformation information from file and
    fits a cubic spline to the data.
"""
from datetime import datetime, time

# External library imports
from scipy import interpolate

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import parsers
from where.lib import config


@plugins.register
def get_gravitational_deformation(rundate):
    """Get excess delay due to gravitational deformation as a function of elevation
    Returns:
        A dictionary of interpolator functions.
    """

    versions = list(config.files.glob_variable("vlbi_gravitational_deformation", "version", r"[\w]+"))

    dates = [datetime.strptime(d, "%Y%b%d") for d in versions]
    if dates:
        max_idx = dates.index(max(dates))
        file_vars = dict(version=versions[max_idx])
    else:
        file_vars = dict()

    parser = parsers.parse_key(file_key="vlbi_gravitational_deformation", file_vars=file_vars)
    data = parser.as_dict() if parser.data_available else dict()

    interpolators = dict()

    for station, values in data.items():
        if (
            datetime.combine(rundate, time.max) > values["start"]
            and datetime.combine(rundate, time.min) < values["end"]
        ):
            interpolators[station] = interpolate.interp1d(values["elevation"], values["delay"], kind="cubic")

    return interpolators
