"""Information about the Lares satellite

Description:
    Useful constants for the Lares satellite.

References:
   https://ilrs.cddis.eosdis.nasa.gov/missions/satellite_missions/current_missions/lars_general.html

"""
# Standard library imports
import math

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.apriori.satellites import satellite


@plugins.register
class Lares(satellite.Slr, satellite.Satellite):

    short_name = "lares"

    # Constants from ILRS webpage [1]:
    cospar_id = 1200601
    area = math.pi * 0.182 ** 2
    mass = 386.8
    radiation_pressure_coefficient = 1.13
