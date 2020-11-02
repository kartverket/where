"""Calculate the delay caused by the antenna axis offset

Description:
------------

Calculate the delay caused by the antenna axis offset.

References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] http://hpiers.obspm.fr/combinaison/documentation/articles/Thermal_Expansion_Modelling_Radio_Telescopes_Nothnagel.pdf

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import log


@plugins.register
def vlbi_axis_offset(dset):
    """Calculate antenna axis offset at both stations

    Args:
        dset (Dataset): Model data.

    Returns:
        Numpy array: Corrections in meters for each observation
    """
    data_out = np.zeros(dset.num_obs)
    for multiplier in dset.for_each_suffix("station"):
        data_out += multiplier * axis_offset_station(dset)

    return data_out


def axis_offset_station(dset):
    """Calculate antenna axis offset at one station

    Args:
        dset:        A Dataset containing model data.

    Returns:
        Numpy array: Antenna axis offset in meters.
    """
    antenna_info = apriori.get("vlbi_antenna_info")
    delays = np.zeros(dset.num_obs)

    sin_a = np.sin(dset.site_pos.azimuth)
    sin_e = np.sin(dset.site_pos.elevation)
    cos_a = np.cos(dset.site_pos.azimuth)
    cos_e = np.cos(dset.site_pos.elevation)
    cos_d = np.cos(dset.src_dir.declination)

    for ivsname in dset.unique("ivsname"):
        site_id = dset.meta[ivsname]["site_id"] if ivsname in dset.meta else ""
        if ivsname not in antenna_info:
            log.warn(f"Missing antenna specifications for ivsname {ivsname!r} ({site_id}). Axis offset correction set to zero.")
            continue

        idx = dset.filter(ivsname=ivsname)
        ao = antenna_info[ivsname]["axis_offset"]
        axis_type = antenna_info[ivsname]["mount"]

        if axis_type == "MO_AZEL":
            delays[idx] = -ao * cos_e[idx]
        elif axis_type == "MO_EQUA":
            delays[idx] = -ao * cos_d[idx]
        elif axis_type == "MO_XYNO":
            delays[idx] = -ao * np.sqrt(1 - (cos_e[idx] * cos_a[idx]) ** 2)
        elif axis_type == "MO_XYEA":
            delays[idx] = -ao * np.sqrt(1 - (cos_e[idx] * sin_a[idx]) ** 2)
        elif axis_type == "MO_RICH":
            # Special case for RICHMOND station
            phi = np.radians(39.06)
            dlam = np.radians(-0.12)
            delays[idx] = -ao * np.sqrt(
                1
                - (
                    sin_e[idx] * np.sin(phi)
                    + cos_e[idx] * np.cos(phi) * (cos_a[idx] * np.cos(dlam) + sin_a[idx] * np.sin(dlam))
                )
                ** 2
            )
        else:
            log.warn(f"Unknown antenna axis type {axis_type!r} for {ivsname} ({site_id}). Correction set to zero")
            continue

    return delays
