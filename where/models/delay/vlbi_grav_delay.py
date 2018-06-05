"""Calculate the gravitational VLBI delay using the consensus model

Description:
------------

Calculate the gravitational delay using the Consensus model as described in the IERS Conventions [1]_, section 11.1.


References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import constant
from where.lib import log
from where.lib import plugins
from where.lib.time import TimeDelta
from where.lib.unit import unit


@plugins.register
def vlbi_grav_delay(dset):
    """Calculate the gravitational delay

    The implementation is described in IERS Conventions [1]_, section 11.1, in particular equation (11.9).

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Gravitational delay in meters for each observation.
    """
    eph = apriori.get("ephemerides", time=dset.time)
    grav_delay = np.zeros(dset.num_obs)

    # List of celestial bodies. Major moons are also recommended, like Titan, Ganymedes, ...
    bodies = [
        "mercury barycenter",
        "venus barycenter",
        "earth",
        "moon",
        "mars barycenter",
        "jupiter barycenter",
        "saturn barycenter",
        "uranus barycenter",
        "neptune barycenter",
        "pluto barycenter",
        "sun",
    ]

    bcrs_vel_earth = eph.vel_bcrs("earth")
    baseline_gcrs = dset.site_pos_2.gcrs_pos - dset.site_pos_1.gcrs_pos
    src_dot_baseline = (dset.src_dir.unit_vector[:, None, :] @ baseline_gcrs[:, :, None])[:, 0, 0]

    # Equation 11.6
    bcrs_site1 = eph.pos_bcrs("earth") + dset.site_pos_1.gcrs_pos
    bcrs_site2 = eph.pos_bcrs("earth") + dset.site_pos_2.gcrs_pos

    for body in bodies:
        try:
            GM_body = constant.get("GM_{}".format(body.split()[0]), source=eph.ephemerides)
        except KeyError:
            log.warn(
                "The GM value of {} is not defined for {}. Correction set to zero.",
                body.split()[0].title(),
                eph.ephemerides,
            )
            continue
        bcrs_body_t1 = eph.pos_bcrs(body)

        # Equation 11.3
        delta_t = TimeDelta(
            np.maximum(0, dset.src_dir.unit_vector[:, None, :] @ (bcrs_body_t1 - bcrs_site1)[:, :, None])[:, 0, 0]
            * unit.second2day
            / constant.c,
            format="jd",
            scale="tdb",
        )
        time_1J = dset.time.tdb - delta_t

        # Equation 11.4
        bcrs_body_t1J = eph.pos_bcrs(body, time=time_1J)
        vector_body_site1 = bcrs_site1 - bcrs_body_t1J

        # Equation 11.5
        vector_body_site2 = bcrs_site2 - bcrs_body_t1J - bcrs_vel_earth / constant.c * src_dot_baseline[:, None]

        # Needed for equation 11.1
        norm_body_site1 = np.linalg.norm(vector_body_site1, axis=1)
        src_dot_vector_body_site1 = (dset.src_dir.unit_vector[:, None, :] @ vector_body_site1[:, :, None])[:, 0, 0]
        nomJ = norm_body_site1 + src_dot_vector_body_site1
        denomJ = np.linalg.norm(vector_body_site2, axis=1) + (
            dset.src_dir.unit_vector[:, None, :] @ vector_body_site2[:, :, None]
        )[
            :, 0, 0
        ]

        # Main correction (equation 11.1)
        grav_delay += 2 * GM_body / constant.c ** 2 * np.log(nomJ / denomJ)

        # Higher order correction  (equation 11.14)
        baseline_dot_vector_body_site1 = (baseline_gcrs[:, None, :] @ vector_body_site1[:, :, None])[:, 0, 0]
        grav_delay += (
            4
            * GM_body
            ** 2
            / constant.c
            ** 4
            * (baseline_dot_vector_body_site1 / norm_body_site1 + src_dot_baseline)
            / (norm_body_site1 + src_dot_vector_body_site1)
            ** 2
        )

    # Denominator (equation 11.9)
    denominator = 1 + (
        (bcrs_vel_earth + dset.site_pos_2.gcrs_vel)[:, None, :] @ dset.src_dir.unit_vector[:, :, None] / constant.c
    )[
        :, 0, 0
    ]

    return grav_delay / denominator
