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
from where.lib import log


@plugins.register
def get_gravitational_deformation(rundate):
    """Get excess delay due to gravitational deformation as a function of elevation
    Returns:
        A dictionary of interpolator functions.
    """

    # Find the file with the highest date in the filename
    versions_1 = list(config.files.glob_variable("vlbi_gravitational_deformation", "version", r"[\w]+"))
    versions_2 = list(config.files.glob_variable("vlbi_gravitational_deformation", "version", r"[\w]+[-][\w]+[-][\w]+"))

    dates_1 = [datetime.strptime(d, "%Y%b%d") for d in versions_1]
    dates_2 = [datetime.strptime(d,"%Y-%m-%d") for d in versions_2]

    # Earlier versions used the date format 2018Nov12
    if dates_1:
        max_date = dates_1.index(max(dates_1))
        file_vars = dict(version=versions_1[max_date])
    else:
        file_vars = dict()

    # Later versions uses the date format 2023-01-20
    if dates_2:
        max_date = dates_2.index(max(dates_2))
        file_vars = dict(version=versions_2[max_date])

    parser = parsers.parse_key(file_key="vlbi_gravitational_deformation", file_vars=file_vars)
    log.debug(f"Using {parser.file_path} as a priori gravitional deformation")
    data = parser.as_dict() if parser.data_available else dict()

    interpolators = dict()

    for station, values in data.items():
        if (
            datetime.combine(rundate, time.max) > values["start"]
            and datetime.combine(rundate, time.min) < values["end"]
        ):
            interpolators[station] = interpolate.interp1d(values["elevation"], values["delay"], kind="cubic")

    return interpolators
