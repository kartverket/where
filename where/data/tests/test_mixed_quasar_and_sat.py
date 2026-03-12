""" Tests for the data.position module"""

# Third party imports
import pytest
import numpy as np

# Where imports
from where.data import time, position, direction

@pytest.fixture
def t():
    v = time.Time(['2024-10-31T18:30:30.000000', '2024-10-31T18:30:30.000000'], scale="utc", fmt="isot")
    return v

@pytest.fixture
def site_pos(t):
    v = position.PosVel(
        [[ 5.08549077e+06,  2.66816176e+06, -2.76869239e+06, -1.70530257e-13,  0.00000000e+00,  5.55111512e-17],
        [-3.94999116e+06,  2.52242125e+06, -4.31170745e+06, 6.39488462e-14,  1.42108547e-14, -1.11022302e-16]],
        system="trs", time=t)
    return v
    
@pytest.fixture
def sat_pos(t):
    v = position.PosVel(
        [[ 1.03030798e+07, -2.07492487e+07, -1.24535559e+07, -3.08877890e+06, -1.41410961e+08,  2.31860459e+08],
        [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]],
        system="trs", time=t)
    return v
    
@pytest.fixture
def src_dir(t):
    v = direction.Direction(
        [[np.nan, np.nan, np.nan],
         [-0.16646546,  0.04565593, -0.98498974]],
        system="gcrs", time=t)
    return v

def test_elevation(site_pos, sat_pos, src_dir):
    e_1_sat = site_pos.elevation_to(sat_pos)[0]
    e_1_src = site_pos.elevation_to(src_dir)[1]
    
    site_pos.other = src_dir
    site_pos.other_2 = sat_pos
    
    e_2_sat = site_pos.elevation[0]
    e_2_src = site_pos.elevation[1]
    
    assert e_1_sat == e_2_sat
    assert e_1_src == e_2_src  

def test_azimuth(site_pos, sat_pos, src_dir):
    a_1_sat = site_pos.azimuth_to(sat_pos)[0]
    a_1_src = site_pos.azimuth_to(src_dir)[1]
    
    site_pos.other = src_dir
    site_pos.other_2 = sat_pos
    
    a_2_sat = site_pos.azimuth[0]
    a_2_src = site_pos.azimuth[1]
    
    assert a_1_sat == a_2_sat
    assert a_1_src == a_2_src

def test_zenith_distance(site_pos, sat_pos, src_dir):
    z_1_sat = site_pos.zenith_distance_to(sat_pos)[0]
    z_1_src = site_pos.zenith_distance_to(src_dir)[1]
    
    site_pos.other = src_dir
    site_pos.other_2 = sat_pos
    
    z_2_sat = site_pos.zenith_distance[0]
    z_2_src = site_pos.zenith_distance[1]
    
    assert z_1_sat == z_2_sat
    assert z_1_src == z_2_src    
    