"""Calculate the gravitational VLBI delay using the consensus model

Description:
------------

Calculate the gravitational delay using the Consensus model as described in the IERS Conventions [1]_, section 11.1.

This model accounts for the gravitational effect of the celestial bodies in the solar system on the delay

References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html





"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori
from midgard.math.constant import constant
from where.lib import log
from where.data.time import TimeDelta


@plugins.register
def vlbi_grav_bodies(dset):
    """Calculate the gravitational delay

    The implementation is described in IERS Conventions [1]_, section 11.1, in particular equation (11.9).

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Gravitational delay in meters for each observation.
    """
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
    output = np.zeros(dset.num_obs)

    # This model is only applicable for far field observations
    idx = ~dset.near_field_obs
    if np.sum(idx) == 0:
        # no far field observations. Skip this model
        return output

    time = dset.time[idx]
    eph = apriori.get("ephemerides", time=time)
    grav_delay = np.zeros(np.sum(idx))
    

    bcrs_vel_earth = eph.vel_bcrs("earth")

    baseline_gcrs = dset.site_pos_2.gcrs.pos[idx] - dset.site_pos_1.gcrs.pos[idx]
    unit_vector_src = dset.src_dir.unit_vector[idx]
    src_dot_baseline = (unit_vector_src[:, None, :] @ baseline_gcrs.mat)[:, 0, 0]
    gcrs_site_vel_2 = dset.site_pos_2.gcrs.vel.val[idx]


    # Equation 11.6
    bcrs_site1 = eph.pos_bcrs("earth") + dset.site_pos_1.gcrs.pos.val[idx]
    bcrs_site2 = eph.pos_bcrs("earth") + dset.site_pos_2.gcrs.pos.val[idx]

    for body in bodies:
        try:
            GM_name = "GM" if body == "earth" else f"GM_{body.split()[0]}"
            GM_body = constant.get(GM_name, source=eph.ephemerides)
        except KeyError:
            log.warn(
                f"The GM value of {body.split()[0].title()} is not defined for {eph.ephemerides}. "
                f"Correction set to zero."
            )
            continue
        bcrs_body_t1 = eph.pos_bcrs(body)

        # Equation 11.3
        delta_t = TimeDelta(
            np.maximum(0, unit_vector_src[:, None, :] @ (bcrs_body_t1 - bcrs_site1)[:, :, None])[:, 0, 0]
            * Unit.second2day
            / constant.c,
            fmt="jd",
            scale="tdb",
        )

        time_1J = time.tdb - delta_t

        # Equation 11.4
        bcrs_body_t1J = eph.pos_bcrs(body, time=time_1J)
        vector_body_site1 = bcrs_site1 - bcrs_body_t1J

        # Equation 11.5
        vector_body_site2 = bcrs_site2 - bcrs_body_t1J - bcrs_vel_earth / constant.c * src_dot_baseline[:, None]

        # Needed for equation 11.1
        norm_body_site1 = np.linalg.norm(vector_body_site1, axis=1)
        src_dot_vector_body_site1 = (unit_vector_src[:, None, :] @ vector_body_site1[:, :, None])[:, 0, 0]
        nomJ = norm_body_site1 + src_dot_vector_body_site1
        denomJ = (
            np.linalg.norm(vector_body_site2, axis=1)
            + (unit_vector_src[:, None, :] @ vector_body_site2[:, :, None])[:, 0, 0]
        )

        # Main correction (equation 11.1)
        grav_delay += 2 * GM_body / constant.c ** 2 * np.log(nomJ / denomJ)

        # Higher order correction  (equation 11.14)
        baseline_dot_vector_body_site1 = (baseline_gcrs.val[:, None, :] @ vector_body_site1[:, :, None])[:, 0, 0]
        grav_delay += (
            4
            * GM_body ** 2
            / constant.c ** 4
            * (baseline_dot_vector_body_site1 / norm_body_site1 + src_dot_baseline)
            / (norm_body_site1 + src_dot_vector_body_site1) ** 2
        )

    # Denominator (equation 11.9)
    denominator = (
        1
        + (
            (bcrs_vel_earth + gcrs_site_vel_2)[:, None, :] @ unit_vector_src[:, :, None]
            / constant.c
        )[:, 0, 0]
    )


    output[idx] = grav_delay / denominator
    return output
