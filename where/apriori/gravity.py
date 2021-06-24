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
import math
from datetime import datetime
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.data import time

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
    GM = gravity_parser.meta["earth_gravity_constant"]
    r = gravity_parser.meta["radius"]

    gr = dict()
    gr["C"] = np.zeros((truncation_level + 1, truncation_level + 1))
    gr["S"] = np.zeros((truncation_level + 1, truncation_level + 1))
    rundate = time.Time(rundate, fmt="datetime", scale="utc")

    if gravity_coefficients["gfc"].shape == ():
        gr["C"][gravity_coefficients["gfc"]["degree"], gravity_coefficients["gfc"]["order"]] = gravity_coefficients[
            "gfc"
        ]["C"]
        gr["S"][gravity_coefficients["gfc"]["degree"], gravity_coefficients["gfc"]["order"]] = gravity_coefficients[
            "gfc"
        ]["S"]
    else:
        for l in gravity_coefficients["gfc"]:
            gr["C"][l["degree"], l["order"]] = l["C"]
            gr["S"][l["degree"], l["order"]] = l["S"]

    if not "gfct" in gravity_coefficients:
        if gravity_field == "EGM2008":
            # Rates not contained in gravity file, apply rates from IERS
            apply_rates(gr, truncation_level, rundate)
        return gr, GM, r

    for l in gravity_coefficients["gfct"]:
        # TODO: Not sure what the xxxx in t0[yyyymmdd.xxxx] actually means (see e.g. EIGEN-6S4)
        # Maybe hours and minutes, but in that case the minutes should have been between 0..59.
        # They are not.
        # It doesn't seem to matter for the orbit determination process,
        # so I set hours = 0 for now

        t0 = time.Time(
            datetime(int(l["t0"][0:4]), int(l["t0"][4:6]), int(l["t0"][6:8]), int(l["t0"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )
        t1 = time.Time(
            datetime(int(l["t1"][0:4]), int(l["t1"][4:6]), int(l["t1"][6:8]), int(l["t1"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )

        if t0 <= rundate < t1:
            gr["C"][l["degree"], l["order"]] += l["C"]
            gr["S"][l["degree"], l["order"]] += l["S"]

    for l in gravity_coefficients["trnd"]:
        t0 = time.Time(
            datetime(int(l["t0"][0:4]), int(l["t0"][4:6]), int(l["t0"][6:8]), int(l["t0"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )
        t1 = time.Time(
            datetime(int(l["t1"][0:4]), int(l["t1"][4:6]), int(l["t1"][6:8]), int(l["t1"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )

        if t0 <= rundate < t1:
            years = (rundate.mjd - t0.mjd) / 365.25
            gr["C"][l["degree"], l["order"]] += l["C"] * years
            gr["S"][l["degree"], l["order"]] += l["S"] * years

    for l in gravity_coefficients["asin"]:
        t0 = time.Time(
            datetime(int(l["t0"][0:4]), int(l["t0"][4:6]), int(l["t0"][6:8]), int(l["t0"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )
        t1 = time.Time(
            datetime(int(l["t1"][0:4]), int(l["t1"][4:6]), int(l["t1"][6:8]), int(l["t1"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )

        if t0 <= rundate < t1:
            years = (rundate.mjd - t0.mjd) / 365.25
            gr["C"][l["degree"], l["order"]] += l["C"] * math.sin(2 * math.pi * years / l["period"])
            gr["S"][l["degree"], l["order"]] += l["S"] * math.sin(2 * math.pi * years / l["period"])

    for l in gravity_coefficients["acos"]:
        t0 = time.Time(
            datetime(int(l["t0"][0:4]), int(l["t0"][4:6]), int(l["t0"][6:8]), int(l["t0"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )
        t1 = time.Time(
            datetime(int(l["t1"][0:4]), int(l["t1"][4:6]), int(l["t1"][6:8]), int(l["t1"][9:11]), 0),
            scale="utc",
            fmt="datetime",
        )

        if t0 <= rundate < t1:
            years = (rundate.mjd - t0.mjd) / 365.25
            gr["C"][l["degree"], l["order"]] += l["C"] * math.cos(2 * math.pi * years / l["period"])
            gr["S"][l["degree"], l["order"]] += l["S"] * math.cos(2 * math.pi * years / l["period"])

    return gr, GM, r


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
    ref_date = time.Time(datetime(2000, 1, 1, 0, 0), fmt="datetime", scale="utc")
    years = (rundate.mjd - ref_date.mjd) / 365.25

    if truncation_level >= 2:
        # Note that the C20- value is the tide free value on page 89 in [4],
        # not the zero tide value in Table 6.2.
        coefficients["C"][2, 0] = -0.484_165_31e-3 + 11.6e-12 * years
    if truncation_level >= 3:
        coefficients["C"][3, 0] += 4.9e-12 * years
    if truncation_level >= 4:
        coefficients["C"][4, 0] += 4.7e-12 * years
