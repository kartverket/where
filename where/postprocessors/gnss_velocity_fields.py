"""Add GNSS site velocity fields to dataset

Description:
------------

Following fields are added to Dataset:

| Field             | Type           | Unit  | Description                                               |
|-------------------|----------------|-------|-----------------------------------------------------------|
| site_vel_east     | numpy.ndarray  |  m/s  | Site velocity east component of topocentric coordinates   |
| site_vel_north    | numpy.ndarray  |  m/s  | Site velocity north component of topocentric coordinates  |
| site_vel_up       | numpy.ndarray  |  m/s  | Site velocity up component of topocentric coordinates     |
| site_vel_h        | numpy.ndarray  |  m/s  | Horizontal site velocity                                  |
| site_vel_v        | numpy.ndarray  |  m/s  | Vertical site velocity                                    |
| site_vel_3d       | numpy.ndarray  |  m/s  | 3D site velocity                                          |
| site_vel_x        | numpy.ndarray  |  m/s  | Site velocity X-coordinate of geocentric coordinates      |
| site_vel_y        | numpy.ndarray  |  m/s  | Site velocity Y-coordinate of geocentric coordinates      |
| site_vel_z        | numpy.ndarray  |  m/s  | Site velocity Z-coordinate of geocentric coordinates      |
| site_vel_sigma_x  | numpy.ndarray  |  m/s  | Standard devication of site velocity X-coordinate         |
| site_vel_sigma_y  | numpy.ndarray  |  m/s  | Standard devication of site velocity Y-coordinate         |
| site_vel_sigma_z  | numpy.ndarray  |  m/s  | Standard devication of site velocity Z-coordinate         |

"""
# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import rotation

# Where imports
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_velocity_fields(dset: "Dataset") -> None:
    """Add GNSS site velocity fields to dataset

    Args:
        dset:     A Dataset containing model data.
    """
    
    # Add velocity in topocentric coordinates
    lat, lon, height = dset.site_pos.pos.llh.T
    vel_enu = np.squeeze(rotation.trs2enu(lat, lon) @  dset.site_vel[:,:,None])
    dset.add_float("site_vel_east", val=vel_enu[:,0], unit="meter/second", write_level="operational")
    dset.add_float("site_vel_north", val=vel_enu[:,1], unit="meter/second", write_level="operational")
    dset.add_float("site_vel_up", val=vel_enu[:,2], unit="meter/second", write_level="operational")
    
    # Add horizontal velocity (HV), vertical velocity (VV) and 3D velocity
    dset.add_float(
            "site_vel_h", 
            val=np.sqrt(vel_enu[:,0] ** 2 + vel_enu[:,1] ** 2),
            unit="meter/second",
            write_level="operational",
    )
    dset.add_float("site_vel_v", val=np.absolute(vel_enu[:,2]), unit="meter/second")
    dset.add_float(
            "site_vel_3d", 
            val=np.sqrt(dset.site_vel[:,0] ** 2 + dset.site_vel[:,1] ** 2 + dset.site_vel[:,2] ** 2),
            #val=np.sqrt(vel_enu[:,0] ** 2 + vel_enu[:,1] ** 2 + vel_enu[:,2] ** 2),
            unit="meter/second",
            write_level="operational",
            
    )

    # Add site velocity and standard deviation of site velocity coordinates
    dset.add_float("site_vel_x", val=dset.site_vel[:, 0], unit="meter/second", write_level="analysis")
    dset.add_float("site_vel_y", val=dset.site_vel[:, 1], unit="meter/second", write_level="analysis")
    dset.add_float("site_vel_z", val=dset.site_vel[:, 2], unit="meter/second", write_level="analysis")
    dset.add_float(
            "site_vel_sigma_x", 
            val=np.sqrt(dset.estimate_cov_site_vel_xx), 
            unit="meter/second",
            write_level="analysis",
    )
    dset.add_float(
            "site_vel_sigma_y", 
            val=np.sqrt(dset.estimate_cov_site_vel_yy), 
            unit="meter/second",
            write_level="analysis",
    )
    dset.add_float(
            "site_vel_sigma_z", 
            val=np.sqrt(dset.estimate_cov_site_vel_zz), 
            unit="meter/second",
            write_level="analysis",
    )

            
