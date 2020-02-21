"""Get gravity coefficients for a given gravity field

Description:

Reads data from gravity files. The gravity field to read is found in the configuration file. Files are assumed to be in
the ICGEM format [1] as available from the International Centre for Global Earth Models (ICGEM) website [2].

Usage of the gravity coefficients is described in the book 'Satellite Orbits' [3] as well as in section 6 of the IERS
Conventions [4]. The IERS Conventions 2010 recommend using the EGM2008 gravity field, which is described in further
detail at [5].

References:
[1] Barthelmes, Franz and Förste, Christoph, The ICGEM-format.
    http://icgem.gfz-potsdam.de/ICGEM/documents/ICGEM-Format-2011.pdf

[2] The International Centre for Global Earth Models (ICGEM).
    http://icgem.gfz-potsdam.de/ICGEM/

[3] Montenbruck, Oliver and Gill, Eberhard, Satellite Orbits,
    Springer Verlag, 2000.

[4] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010).
    http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

[5] http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/

"""

# Standard library imports
from datetime import datetime

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log
from where import parsers


@plugins.register
def get_gravity_coefficients(gravity_field, truncation_level, rundate):
    """Get coefficients for a gravity field

    The coefficient files come with normalized gravity coefficients.

    The coefficients are returned in a dict with keys `C` and `S`,
    each containing a table of coefficients `C[i, j]` and `S[i, j]`.

    Args:
        gravity_field:    Name of gravity field.
        truncation_level: Level of degree and order where coeffs are truncated.
        rundate:          Start time of integration.

    Returns:
        A dictionary containing C and S coefficients.
    """
    log.info(f"Reading gravity field {gravity_field} up to degree and order {truncation_level}")
    try:
        gravity_parser = parsers.parse_key(
            "gravity_coefficients", file_vars=dict(gravity_field=gravity_field), num_degrees=truncation_level
        )
    except FileNotFoundError:
        log.fatal(f"Unknown gravity field {gravity_field}, exiting!")

    gravity_coefficients = gravity_parser.as_dict()
    if gravity_field == "egm2008":  # TODO: Why do we treat this differently ... Should be done by properties not name
        apply_rates(gravity_coefficients, truncation_level, rundate)

    return gravity_coefficients


def apply_rates(coefficients, truncation_level, rundate):
    """Apply rates to coefficients for the gravitational field

    See IERS Conventions [4], equation (6.4).

    Args:
        coefficients:     Dict containing tables of coefficients.
        truncation_level: Level of degree and order where coeffs are truncated.
        rundate:          Start time of integration.

    Returns:
        None. Coefficients are changed in place.
    """

    years = (rundate - datetime(2000, 1, 1)).days / 365.25

    if truncation_level >= 2:
        # Note that the C20- value is the tide free value on page 89 in [4],
        # not the zero tide value in Table 6.2.
        coefficients["C"][2, 0] = -0.484_165_31e-3 + 11.6e-12 * years
    if truncation_level >= 3:
        coefficients["C"][3, 0] += 4.9e-12 * years
    if truncation_level >= 4:
        coefficients["C"][4, 0] += 4.7e-12 * years
