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
def src_dir(t):
    v = direction.Direction(
        [[np.nan, np.nan, np.nan],
         [-0.16646546,  0.04565593, -0.98498974]],
        system="gcrs", time=t)
    return v

@pytest.fixture
def sat_pos(t):
    v = position.PosVel(
        [[ 1.03030798e+07, -2.07492487e+07, -1.24535559e+07, -3.08877890e+06, -1.41410961e+08,  2.31860459e+08],
        [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]],
        system="trs", time=t)
    return v

@pytest.fixture
def sat_pos_kepler(t):
    v = position.PosVel(
        [[ 25015181.018465, 0.707977170873199, 0.121662175957290, 3.024483909022929, 1.597899323919624, 2.772570719534964],
        [ 25015182.018465, 0.717977170873199, 0.131662175957290, 3.124483909022929, 1.697899323919624, 2.872570719534964]],
        system="kepler", time=t)
    return v

@pytest.fixture
def site_pos(t, src_dir, sat_pos):
    v = position.PosVel(
        [[ 5.08549077e+06,  2.66816176e+06, -2.76869239e+06, -1.70530257e-13,  0.00000000e+00,  5.55111512e-17],
        [-3.94999116e+06,  2.52242125e+06, -4.31170745e+06, 6.39488462e-14,  1.42108547e-14, -1.11022302e-16]],
        system="trs", time=t, other=src_dir, other_2=sat_pos)
    return v

@pytest.fixture
def site_pos_gcrs(t, src_dir, sat_pos):
    v = position.PosVel(
        [[ 5.08549077e+06,  2.66816176e+06, -2.76869239e+06, -1.70530257e-13,  0.00000000e+00,  5.55111512e-17],
        [-3.94999116e+06,  2.52242125e+06, -4.31170745e+06, 6.39488462e-14,  1.42108547e-14, -1.11022302e-16]],
        system="gcrs", time=t, other=src_dir, other_2=sat_pos)
    return v

@pytest.fixture
def delta_pos(t, site_pos):
    v = position.PosVelDelta(
        [[ 0.1,  0.2, 0.3, 0.01, 0.02, 0.03],
        [0.11,  0.22, 0.33, 0.011,  0.022, 0.033]],
        system="trs", time=t, ref_pos=site_pos)
    return v

@pytest.fixture
def delta_pos_yaw(t, sat_pos):
    v = position.PosVelDelta(
        [[ 0.1,  0.2, 0.3, 0.01, 0.02, 0.03],
        [0.11,  0.22, 0.33, 0.011,  0.022, 0.033]],
        system="yaw", time=t, ref_pos=sat_pos)
    return v
    
