"""Information about the Etalon satellites

Description:
    Useful constants for the Etalon 1 and 2 satellites.

References:
   [1] https://ilrs.cddis.eosdis.nasa.gov/missions/satellite_missions/current_missions/eta1_general.html
   [2] https://ilrs.cddis.eosdis.nasa.gov/missions/satellite_missions/current_missions/eta2_general.html
   [3] M.H.Torrence, P.J.Dunn, R.Kolenkiewicz: Characteristics of the LAGEOS and ETALON satellites orbits,
       Advances of Space Research (1995), Volume 16, Issue 12, pp 21 - 24.

"""
# Standard library imports
import math

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.apriori.satellites import satellite


@plugins.register
class Etalon1(satellite.Slr, satellite.Satellite):

    short_name = "eta1"

    # Constants from ILRS webpage [1]:
    cospar_id = 8900103
    area = math.pi * (1.294 / 2) ** 2
    mass = 1415
    radiation_pressure_coefficient = 1.22


@plugins.register
class Etalon2(satellite.Slr, satellite.Satellite):

    short_name = "eta2"

    # Constants from ILRS webpage [1]:
    cospar_id = 8903903
    area = math.pi * (1.294 / 2) ** 2
    mass = 1415
    radiation_pressure_coefficient = 1.25
