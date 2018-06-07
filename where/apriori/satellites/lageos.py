"""Information about the Lageos satellites

Description:
    Useful constants for the Lageos 1 and 2 satellites.

References:
   [1] http://ilrs.gsfc.nasa.gov/missions/satellite_missions/current_missions/lag1_general.html

   [2] D.P.Rubincam: On the Secular Decrease in Semimajor Axis of Lageos's Orbit, Celestial Mechanics 26 (1982),
             pp 361-382.



"""
# Standard library imports
import math

# Where imports
from where.lib import plugins
from where.apriori.satellites import satellite


@plugins.register
class Lageos1(satellite.Slr, satellite.Satellite):

    short_name = "lag1"

    # Constants from ILRS webpage [1]:
    cospar_id = 7603901
    area = math.pi * 0.30 ** 2
    mass = 406.965
    radiation_pressure_coefficient = 1.13

    # @todo Read these values from apriori file in the future?
    # Easier to estimate drag coefficient times atmospheric density.
    # Initial values: constants from Rubincam [1].
    rho = 1.5e-18  # rho on page 376
    drag_coefficient = 2.2  # C_D on page 376
    drag_product = rho * drag_coefficient


@plugins.register
class Lageos2(satellite.Slr, satellite.Satellite):

    short_name = "lag2"

    # Constants from ILRS webpage [1]:
    cospar_id = 9207002
    area = math.pi * 0.30 ** 2
    mass = 405.38
    radiation_pressure_coefficient = 1.13

    # @todo Read these values from apriori file in the future?
    # Easier to estimate drag coefficient times atmospheric density.
    # Initial values: constants from Rubincam [1].
    rho = 1.5e-18  # rho on page 376
    drag_coefficient = 2.2  # C_D on page 376
    drag_product = rho * drag_coefficient
