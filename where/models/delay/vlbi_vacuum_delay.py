"""Calculate the VLBI geometric vacuum delay using the consensus model

Description:
------------

Calculate the geometric vacuum delay using the Consensus model as described in the IERS Conventions :cite:`iers2010`,
section 11.1. This takes into account both the rotation and the velocity of the Earth, as well as the gravitational
potential. Apriori Earth Orientation Parameters are interpolated and corrected for subdaily tides.





"""

# External library imports
import numpy as np

# Where imports
from where import apriori
from midgard.math.constant import constant
from where.lib import plugins


@plugins.register
def vlbi_vacuum_delay(dset):
    r"""Calculate the theoretical delay dependent on the baseline

    The implementation is described in IERS Conventions :cite:`iers2010`, section 11.1, in particular equation
    (11.9). We do not take the gravitational delay into account here (see
    :mod:`where.models.delay.vlbi_gravitational_delay`), and multiply by :math:`c` to get the correction in
    meters. Thus, we implement the following equation:

    .. math::
       \mathrm{correction} = \frac{- \hat K \cdot \vec b \bigl[ 1 - \frac{(1 + \gamma) U}{c^2}
                             - \frac{| \vec V_\oplus |^2}{2 c^2} - \frac{\vec V_\oplus \cdot \vec w_2}{c^2} \bigr]
                             - \frac{\vec V_\oplus \cdot \vec b}{c} \bigl[ 1
                             + \frac{\hat K \cdot \vec V_\oplus}{2 c} \bigr]}{1
                             + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}}

    with

    * :math:`\hat K` -- the unit vector from the barycenter to the source in the absence of gravitational or
      aberrational bending,

    * :math:`\vec b` -- the GCRS baseline vector at the time :math:`t_1` of arrival, :math:`\vec x_2(t_1) - \vec
      x_1(t_1)`,

    * :math:`\gamma` -- the parameterized post-Newtonian (PPN) gamma, equal to 1 in general relativity theory,

    * :math:`U` -- the gravitational potential at the geocenter, neglecting the effects of the Earth's mass. At the
      picosecond level, only the solar potential need be in included in :math:`U` so that :math:`U = G M_\odot / | \vec
      R_{\oplus_\odot} |` where :math:`\vec R_{\oplus_\odot}` is the vector from the Sun to the geocenter,

    * :math:`\vec V_\oplus` -- the barycentric velocity of the geocenter,

    * :math:`\vec w_2` -- the geocentric velocity of station 2.

    Each term in the correction is calculated in separate functions. and stored in the Dataset in a table called
    ``vlbi_vacuum_delay``.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array: Projected baseline in meters for each observation.

    """
    eph = apriori.get("ephemerides", time=dset.time)
    terms = (term_1, term_2, term_3, term_4, term_5, term_6)

    # Calculate each term of the geometric vacuum delay
    vel_earth = eph.vel_bcrs("earth")
    baseline_gcrs = dset.site_pos_2.gcrs_pos - dset.site_pos_1.gcrs_pos
    denominator = (
        1
        + ((vel_earth + dset.site_pos_2.gcrs_vel)[:, None, :] @ dset.src_dir.unit_vector[:, :, None] / constant.c)[
            :, 0, 0
        ]
    )
    proj_Kb = (dset.src_dir.unit_vector[:, None, :] @ baseline_gcrs[:, :, None])[:, 0, 0] / denominator
    proj_Vb = (vel_earth[:, None, :] @ baseline_gcrs[:, :, None] / constant.c)[:, 0, 0] / denominator

    for term_func in terms:
        field = "vlbi_vacuum_delay_" + term_func.__name__
        values = term_func(dset, proj_Kb, proj_Vb, vel_earth)
        if field in dset.fields:
            dset[field][:] = values
        else:
            dset.add_float(field, table="vlbi_vacuum_delay", val=values, write_level="detail")

    return np.sum(dset.get_table("vlbi_vacuum_delay"), axis=1)


def term_1(dset, proj_Kb, _, _ve):
    r"""Main part of the vacuum delay is the baseline in the source direction

    The term :math:`\hat K \cdot \vec b \cdot \bigl( -1 \bigr)` scaled by the denominator :math:`1 + \frac{\hat K \cdot
    (\vec V_\oplus + \vec w_2)}{c}`.

    Args:
        dset:    Model input data.
        proj_Kb: Scaled projection of baseline in direction of source.

    Returns:
        Numpy array: Part of vacuum delay.
    """
    return -proj_Kb


def term_2(dset, proj_Kb, _, _ve):
    r"""Part of the vacuum delay dependent on the gravitational potential

    The term :math:`\hat K \cdot \vec b \cdot \frac{(1 + \gamma) U}{c^2}` scaled by the denominator :math:`1 +
    \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.

    The parameterized post-Newtonian (PPN) gamma, :math:`\gamma` is equal to 1 in general relativity theory.

    The gravitational potential at the geocenter, \f$ U \f$, neglecting the effects of the Earth's mass is
    calculated. Following table 11.1 in IERS Conventions [2], only the solar potential need to be included at the
    picosecond level. That is

    \f[ U = G M_\odot / | \vec R_{\oplus_\odot} | \f]

    where \f$ \vec R_{\oplus_\odot} \f$ is the vector from the Sun to the geocenter. We calculate the latter using the
    ephemerides.

    Args:
        dset:    Model input data.
        proj_Kb: Scaled projection of baseline in direction of source.

    Returns:
        Numpy array: Part of vacuum delay.
    """
    gamma = 1.0
    eph = apriori.get("ephemerides", time=dset.time)
    grav_potential = constant.GM_sun / np.linalg.norm(eph.pos_gcrs("sun"), axis=1)
    return proj_Kb * (1 + gamma) * grav_potential / constant.c ** 2


def term_3(dset, proj_Kb, _, vel_earth):
    r"""Correction to delay based on earth's movement in space

    The term :math:`\hat K \cdot \vec b \cdot \frac{| \vec V_\oplus |^2}{2 c^2}` scaled by the denominator :math:`1 +
    \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.

    Args:
        dset:    Model input data.
        proj_Kb: Scaled projection of baseline in direction of source.

    Returns:
        Numpy array: Part of vacuum delay.
    """
    return proj_Kb * 0.5 * (vel_earth[:, None, :] @ vel_earth[:, :, None] / constant.c ** 2)[:, 0, 0]


def term_4(dset, proj_Kb, _, vel_earth):
    r"""Correction to the delay caused by earth's rotation

    The term :math:`\hat K \cdot \vec b \cdot \frac{\vec V_\oplus \cdot \vec w_2}{c^2}` scaled by the denominator
    :math:`1 + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.

    Args:
        dset:    Model input data.
        proj_Kb: Scaled projection of baseline in direction of source.

    Returns:
        Numpy array: Part of vacuum delay.
    """
    return proj_Kb * (vel_earth[:, None, :] @ dset.site_pos_2.gcrs_vel[:, :, None] / constant.c ** 2)[:, 0, 0]


def term_5(dset, _, proj_Vb, _ve):
    r"""Part of the delay due to earth's movement in space

    The term :math:`- \frac{\vec V_\oplus \cdot \vec b}{c} \cdot \bigl( 1 \bigr)` scaled by the denominator :math:`1 +
    \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.

    Args:
        dset:    Model input data.
        proj_Vb: Scaled projection of baseline in direction of earth's movement.

    Returns:
        Numpy array: Part of vacuum delay.
    """
    return -proj_Vb


def term_6(dset, _, proj_Vb, vel_earth):
    r"""Correction to earth's movement in space

    The term :math:`- \frac{\vec V_\oplus \cdot \vec b}{c} \cdot \frac{\hat K \cdot \vec V_\oplus}{2 c}` scaled by the
    denominator :math:`1 + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}`.

    Args:
        dset:    Model input data.
        proj_Vb: Scaled projection of baseline in direction of earth's movement.

    Returns:
        Numpy array: Part of vacuum delay.
    """
    return -proj_Vb * 0.5 * (dset.src_dir.unit_vector[:, None, :] @ vel_earth[:, :, None] / constant.c)[:, 0, 0]
