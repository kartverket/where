"""Correct for effect of GNSS carrier phase wind-up


TODO: replace get_yaw_coord_sys() by functions given in posvel_table.py

$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External library imports
import numpy as np

# Where imports
from where import apriori
from where.lib import plugins


@plugins.register
def gnss_carrier_phase_wind_up(dset):
    """Determine carrier phase wind-up correction

    The correction is caluclated after the description in Section 5.5 in :cite:`subirana2013` based on :cite:`wu1993`.

    Args:
        dset (Dataset):   Model data.

    Returns:
        numpy.ndarray:    Carrier phase wind-up correction
    """
    correction = np.zeros(dset.num_obs)

    # Loop over all observed satellites
    for sat in dset.unique("satellite"):

        idx = dset.filter(satellite=sat)
        num_sat_obs = len(dset.sat_posvel.itrs[idx])
        phi = np.zeros(num_sat_obs)

        # Get orientation of satellite transmitter antenna, which are unit vectors from yaw-steering satellite
        # coordinate system
        unit_x, unit_y, unit_z = get_yaw_coord_sys(dset.time[idx], dset.sat_posvel.itrs_pos[idx])

        # Get orientation of receiver by defining East, North and Up unit vectors
        unit_e = dset.site_pos.enu_east[idx]
        unit_n = dset.site_pos.enu_north[idx]
        unit_u = dset.site_pos.enu_up[idx]  # TODO: This variable is never used ...

        # Determine unit line-of-sight vector pointing from transmitter to the receiver
        sat_pos = -dset.sat_posvel.itrs_pos[idx]
        sat_pos_norm = np.linalg.norm(sat_pos, axis=1)
        unit_sat_rec = sat_pos / sat_pos_norm[:, None]

        # Determine effective dipole of receiver and satellite antenna
        #
        # Size of matrices:
        # (nobs,3) = (nobs,3) - (nobs,3) * ((nobs,1,3) * (nobs,3,1)) + (nobs,3) x (nobs,3)
        dipole_rec = (
            unit_e
            - unit_sat_rec
            * (unit_sat_rec[:, None, :] @ unit_e[:, :, None])[:, :, 0]
            + np.cross(unit_sat_rec, unit_n)
        )
        dipole_sat = (
            unit_x
            + unit_sat_rec
            * (unit_sat_rec[:, None, :] @ unit_x[:, :, None])[:, :, 0]
            - np.cross(unit_sat_rec, unit_y)
        )
        dipole_rec_norm = np.linalg.norm(dipole_rec, axis=1)
        dipole_sat_norm = np.linalg.norm(dipole_sat, axis=1)

        zeta = (unit_sat_rec[:, None, :] @ np.cross(dipole_sat, dipole_rec)[:, :, None])[:, :, 0]
        dphi = np.sign(zeta) * np.arccos(
            (dipole_sat[:, None, :] @ dipole_rec[:, :, None])[:, :, 0] / (dipole_sat_norm * dipole_rec_norm)[:, None]
        )

        phi[0] = dphi[0, 0]
        for ii in range(1, num_sat_obs):

            # NOTE to np.rint(): For values exactly halfway between rounded decimal values, Numpy rounds to the nearest
            #      even value. Thus 1.5 and 2.5 round to 2.0, -0.5 and 0.5 round to 0.0, etc..
            phi[ii] = dphi[ii] + 2 * np.pi * np.rint((phi[ii - 1] - dphi[ii]) / (2 * np.pi))

        correction[idx] = phi

    return -correction


def get_yaw_coord_sys(time, sat_pos):
    """

    Args:
        time (where.lib.time.Time):   Where Time object.
        sat_pos (numpy.ndarray):      Satellite position in ITRS.

    Returns:
        tuple:  with following `numpy.ndarray` unit vectors defining the yaw-steering coordinate axes

        ===============  ============================================================================================
         Elements         Description
        ===============  ============================================================================================
         unit_x           Unit vector of x-axis lying in the Earth-Satellite-Sun plane in [m]
         unit_y           Unit vector of y-axis, which is the normal vector of the Earth-Satellite-Sun plane in [m]
         unit_z           Unit vector of z-axis pointing to the Earth's center
        ===============  ============================================================================================

    """

    num_obs = len(time)
    unit_y = np.zeros((num_obs, 3))
    unit_z = np.zeros((num_obs, 3))

    sat_pos_norm = np.linalg.norm(sat_pos, axis=1)
    unit_z = -sat_pos / sat_pos_norm[:, None]

    unit_sat_sun_pos = get_satellite_sun_vector(time, sat_pos)
    y = np.cross(unit_sat_sun_pos, sat_pos)
    y_norm = np.linalg.norm(y, axis=1)
    unit_y = y / y_norm[:, None]

    unit_x = np.cross(unit_y, unit_z)

    return unit_x, unit_y, unit_z


def get_satellite_sun_vector(time, sat_pos):
    """Determine unit vector pointing from satellite to Sun in ITRS

    The determination of the vector pointing from satellite to Sun is based on Eq. 5.77 in :cite:`subirana2013`.

    Args:
        time (where.lib.time.Time):     Where Time object.
        sat_pos (numpy.ndarray):        Satellite position in ITRS.

    Returns:
        numpy.ndarray:  unit vectors pointing from satellite to Sun in ITRS and in unit of meters
    """

    num_obs = len(time)
    unit_sat_sun_pos = np.zeros((num_obs, 3))

    # Get Sun position vector in ITRS
    eph = apriori.get("ephemerides", time=time)

    # TODO:
    # Actually the JPL ephemeris are given in the BCRS with Solar System barycenter as origin and not Earth mass
    # center.  So in principle the sun position vector has to be transformed from the BCRS to the GCRS. What are the
    # consequences, if we do not consider these corrections?
    sun_pos_itrs = eph.pos_itrs("sun")

    # Determination of vector between satellite and Sun
    sat_sun_pos = sun_pos_itrs - sat_pos
    sat_sun_pos_norm = np.linalg.norm(sat_sun_pos, axis=1)
    unit_sat_sun_pos = sat_sun_pos / sat_sun_pos_norm[:, None]

    return unit_sat_sun_pos
