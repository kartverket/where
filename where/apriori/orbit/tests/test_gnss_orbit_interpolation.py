""" Test GNSS orbit interpolation methods.

Description:
------------
For the testing of GNSS orbit interpolation a precise orbit dataset with 15 min data interval is compared against a
precise orbit dataset with 5 min data interval. The comparison follows the example described in Schenewerk
:cite:'schenewerk2002' and the test orbit data are downloaded from `http://www.ngs.noaa.gov/gps-toolbox/sp3intrp.htm`.

-------

$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""


# Standard library imports
from datetime import datetime

# External library imports
import numpy as np
import matplotlib.pyplot as plt
import pytest

# Where imports
from where import parsers
from where import data
from where.lib import config
from where.lib import files
from where.lib import mathp


def _plot(x_val, y_val, title=""):
    """Plot

    Args:
        x_val       - Array with x-values
        y_val       - Array with y-values
        title       - title
    """

    plt.title(title, fontsize=12)
    plt.plot(x_val, y_val, "bo")
    plt.xlabel("Time")
    plt.legend(fontsize=8, loc=1)
    plt.show()


# Interpolation settings
moving_window = True
interpolation_method = "lagrange"  # Methods: 'lagrange', 'interp1d' or 'InterpolatedUnivariateSpline'
window_size = 10

# Definition of dummy date
year = 2002
month = 1
day = 1
hour = 0
minute = 0
second = 0
rundate = datetime(year, month, day, hour, minute, second)
file_vars = config.date_vars(rundate)

# Define 15 min dataset
file_path = files.path(file_key="test_gnss_orbit_interpolation_15min", file_vars=file_vars)
orb_15min = data.Dataset(
    rundate, tech=None, stage=None, dataset_name="gnss_precise_orbit_15min", dataset_id=0, empty=True
)
parser = parsers.parse(parser_name="orbit_sp3c", file_path=file_path, rundate=rundate)
parser.write_to_dataset(orb_15min)

# Define 5 min control dataset
file_path = files.path(file_key="test_gnss_orbit_interpolation_5min", file_vars=file_vars)
orb_5min = data.Dataset(
    rundate, tech=None, stage=None, dataset_name="gnss_precise_orbit_15min", dataset_id=0, empty=True
)
parser = parsers.parse(parser_name="orbit_sp3c", file_path=file_path, rundate=rundate)
parser.write_to_dataset(orb_5min)

# Interpolation timestamps
obs_time = orb_5min.time.gps.gpssec
sat_pos = np.zeros((orb_5min.num_obs, 3))
sat_vel = np.zeros((orb_5min.num_obs, 3))


@pytest.mark.skip(reason="Test not fully implemented yet")
def test_gnss_orbit_interpolation():
    for sat in orb_5min.unique("satellite"):

        if sat == "G01":

            orb_idx = orb_15min.filter(satellite=sat)
            obs_idx = orb_5min.filter(satellite=sat)

            if moving_window:
                sat_pos[obs_idx], sat_vel[obs_idx] = mathp.moving_window_interpolation(
                    orb_15min.time.gps.gpssec[orb_idx],
                    orb_15min.sat_pos.itrs[orb_idx],
                    orb_5min.time.gps.gpssec[obs_idx],
                    model=interpolation_method,
                    window_size=window_size,
                )
            else:
                sat_pos[obs_idx], sat_vel[obs_idx] = mathp.interpolation(
                    orb_15min.time.gps.gpssec[orb_idx],
                    orb_15min.sat_pos.itrs[orb_idx],
                    orb_5min.time.gps.gpssec[obs_idx],
                    model=interpolation_method,
                )

            diff_orb = sat_pos[obs_idx] - orb_5min.sat_pos.itrs[obs_idx]
            diff_orb_3d = np.linalg.norm(diff_orb, axis=1)
            rms_3d = np.sqrt(np.mean(np.square(diff_orb_3d)))

            _plot(
                orb_5min.time.gps.gpssec[obs_idx],
                diff_orb_3d,
                title="3D-difference between interpolated and control orbit values\n for satellite {:s} (RMS: {:>6.4f} m)".format(
                    sat, rms_3d
                ),
            )
