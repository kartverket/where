"""Adds dilution of precision (DOP) to dataset

Description:
------------
Dilution of precision calculation is based on estimated covariance matrix of unknowns. GDOP, PDOP, HDOP, VDOP and TDOP
is added to dataset.

TODO: Check if the calculation of HDOP and VDOP is correct. 
"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_dop_cov(dset: "Dataset") -> None:
    """Adds dilution of precision (DOP) to dataset

    Args:
        dset:     A Dataset containing model data.
    """

    # PDOP
    pdop = (
        np.sqrt(dset.estimate_cov_site_pos_xx + dset.estimate_cov_site_pos_yy + dset.estimate_cov_site_pos_zz)
        / dset.estimate_variance_factor
    )

    if "pdop" in dset.fields:
        dset["pdop"][:] = pdop
        log.debug(f"{_SECTION}: Update pdop field in Dataset.")

    else:
        dset.add_float("pdop", val=pdop)

    # GDOP
    gdop = np.sqrt(
        dset.estimate_cov_site_pos_xx
        + dset.estimate_cov_site_pos_yy
        + dset.estimate_cov_site_pos_zz
        + dset.estimate_cov_rcv_clock_tt
    ) / (dset.estimate_variance_factor)

    if "gdop" in dset.fields:
        dset["gdop"][:] = pdop
        log.debug(f"{_SECTION}: Update gdop field in Dataset.")

    else:
        dset.add_float("gdop", val=gdop)

    # TDOP
    tdop = np.sqrt(dset.estimate_cov_rcv_clock_tt) / dset.estimate_variance_factor

    if "tdop" in dset.fields:
        dset["tdop"][:] = pdop
        log.debug(f"{_SECTION}: Update tdop field in Dataset.")

    else:
        dset.add_float("tdop", val=tdop)

    # HDOP and VDOP
    #
    # Epochwise estimation or over whole time period
    dop_xyz = True
    if dop_xyz:

        # HDOP (xyz)
        hdop = np.sqrt(dset.estimate_cov_site_pos_xx + dset.estimate_cov_site_pos_yy) / dset.estimate_variance_factor

        # VDOP
        vdop = np.sqrt(dset.estimate_cov_site_pos_zz) / dset.estimate_variance_factor

    else:
        if config.tech.estimate_epochwise.bool:
            hdop = np.zeros(dset.num_obs)
            vdop = np.zeros(dset.num_obs)

            for epoch in sorted(set(dset.time.gps.mjd)):
                cov_xyz = np.zeros((3, 3))
                idx = dset.time.gps.mjd == epoch

                cov_xyz[0][0] = dset.estimate_cov_site_pos_xx[idx][0]
                cov_xyz[0][1] = dset.estimate_cov_site_pos_xy[idx][0]
                cov_xyz[1][0] = dset.estimate_cov_site_pos_xy[idx][0]
                cov_xyz[0][2] = dset.estimate_cov_site_pos_xz[idx][0]
                cov_xyz[2][0] = dset.estimate_cov_site_pos_xz[idx][0]
                cov_xyz[1][1] = dset.estimate_cov_site_pos_yy[idx][0]
                cov_xyz[1][2] = dset.estimate_cov_site_pos_yz[idx][0]
                cov_xyz[2][1] = dset.estimate_cov_site_pos_yz[idx][0]
                cov_xyz[2][2] = dset.estimate_cov_site_pos_zz[idx][0]

                R = dset.site_pos._enu2itrs[idx][0]
                sigma0 = dset.estimate_variance_factor[idx][0]
                q_enu = R.T @ (cov_xyz / sigma0) @ R
                hdop[idx] = np.sqrt(q_enu[0][0] + q_enu[1][1])
                vdop[idx] = np.sqrt(q_enu[2][2])

        else:
            cov_xyz = np.zeros((3, 3))
            cov_xyz[0][0] = dset.estimate_cov_site_pos_xx[0]
            cov_xyz[0][1] = dset.estimate_cov_site_pos_xy[0]
            cov_xyz[1][0] = dset.estimate_cov_site_pos_xy[0]
            cov_xyz[0][2] = dset.estimate_cov_site_pos_xz[0]
            cov_xyz[2][0] = dset.estimate_cov_site_pos_xz[0]
            cov_xyz[1][1] = dset.estimate_cov_site_pos_yy[0]
            cov_xyz[1][2] = dset.estimate_cov_site_pos_yz[0]
            cov_xyz[2][1] = dset.estimate_cov_site_pos_yz[0]
            cov_xyz[2][2] = dset.estimate_cov_site_pos_zz[0]

            R = dset.site_pos._enu2itrs[0]
            sigma0 = dset.estimate_variance_factor[0]
            q_enu = R.T @ (cov_xyz / sigma0) @ R
            hdop = np.repeat(np.sqrt(q_enu[0][0] + q_enu[1][1]), dset.num_obs)
            vdop = np.repeat(np.sqrt(q_enu[2][2]), dset.num_obs)

    if "hdop" in dset.fields:
        dset["hdop"][:] = hdop
        log.debug(f"{_SECTION}: Update hdop field in Dataset.")

    else:
        dset.add_float("hdop", val=hdop)

    if "vdop" in dset.fields:
        dset["vdop"][:] = vdop
        log.debug(f"{_SECTION}: Update vdop field in Dataset.")

    else:
        dset.add_float("vdop", val=vdop)
    dset.add_float
