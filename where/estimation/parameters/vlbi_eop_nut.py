"""Calculate the partial derivatives of the celestial pole offset Earth Orientation Parameters.

Description:
------------

Calculate the partial derivatives of the :math:`X` and :math:`Y` Earth orientation parameters.

This is done according to equations (2.37) - (2.46) in Teke :cite:`teke2011`.



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import plugins
from where.lib import rotation
from where.ext import sofa_wrapper as sofa
from where.lib.unit import unit


@plugins.register
def eop_nut(dset):
    """Calculate the partial derivative of the celestial pole offset Earth Orientation Parameters

    Args:
        data:     A Dataset containing model data.

    Returns:
        Tuple: Array of partial derivatives, and list of names of derivatives
    """
    column_names = ["x", "y"]
    src_dir = dset.src_dir.unit_vector[:, None, :]
    baseline = (dset.site_pos_2.itrs_pos - dset.site_pos_1.itrs_pos)[:, :, None]

    # Read EOP values
    eop = apriori.get("eop", time=dset.time)

    # CIP and CIO
    X, Y = sofa.vectorized_xy06(dset.time)
    s = sofa.vectorized_s06(dset.time, X, Y)

    # Add Celestial Intermediate Pole corrections
    X += eop.dx * unit.arcsec2rad
    Y += eop.dy * unit.arcsec2rad

    partials = np.zeros((dset.num_obs, 2))

    E = np.arctan2(Y, X)
    d = np.arccos(np.sqrt(1 - (X ** 2 + Y ** 2)))

    # Some intermediate partials
    dE_dX = (-Y / (X ** 2 + Y ** 2))[:, None, None]
    dmE_dX = (Y / (X ** 2 + Y ** 2))[:, None, None]
    ds_dX = (-Y / 2)[:, None, None]
    dmd_dX = (X / (np.sqrt(1 - (X ** 2 + Y ** 2)) * np.sqrt(X ** 2 + Y ** 2)))[:, None, None]

    dE_dY = (X / (X ** 2 + Y ** 2))[:, None, None]
    dmE_dY = (-X / (X ** 2 + Y ** 2))[:, None, None]
    ds_dY = (-X / 2)[:, None, None]
    dmd_dY = (Y / (np.sqrt(1 + (X ** 2 + Y ** 2)) * np.sqrt(X ** 2 + Y ** 2)))[:, None, None]

    # Rotation matrices
    R3_E = rotation.R3(E)
    R3_s = rotation.R3(s)
    R2_md = rotation.R2(-d)
    R3_mE = rotation.R3(-E)
    dR3_s = rotation.dR3(s)
    dR3_E = rotation.dR3(E)
    dR3_mE = rotation.dR3(-E)
    dR2_md = rotation.dR2(-d)

    # Celestial pole offset X
    dQ_dX = (
        dmE_dX
        * dR3_mE @ R2_md @ R3_E @ R3_s
        + dmd_dX
        * R3_mE @ dR2_md @ R3_E @ R3_s
        + dE_dX
        * R3_mE @ R2_md @ dR3_E @ R3_s
        + ds_dX
        * R3_mE @ R2_md @ R3_E @ dR3_s
    )
    partials[:, 0] = (src_dir @ dQ_dX @ sofa.R(dset.time) @ sofa.W(dset.time) @ baseline)[:, 0, 0]

    # Celestial pole offset Y
    dQ_dY = (
        dmE_dY
        * dR3_mE @ R2_md @ R3_E @ R3_s
        + dmd_dY
        * R3_mE @ dR2_md @ R3_E @ R3_s
        + dE_dY
        * R3_mE @ R2_md @ dR3_E @ R3_s
        + ds_dY
        * R3_mE @ R2_md @ R3_E @ dR3_s
    )
    partials[:, 1] = (src_dir @ dQ_dY @ sofa.R(dset.time) @ sofa.W(dset.time) @ baseline)[:, 0, 0]

    return partials, column_names, "meter per radian"
