# TODO move to test_position.py? create test_directon.py and and test_time.py?
# similar tests for pos, posdelta, posveldelta?
import pytest

# TODO can this import be avoided?
from conftest import *

@pytest.fixture()
def p(request):
    return request.getfixturevalue(request.param.__name__)

@pytest.mark.parametrize("p", (site_pos, delta_pos, delta_pos_yaw, site_pos_gcrs, sat_pos_kepler, src_dir), indirect=True)
def test_sliced_fields(p):

    sliced_p_1 = p[0] # Pick one element
    sliced_p_2 = p[0:1] # Pick one slice (of size 1)
    for f in p.fieldnames():
        field = getattr(p, f)
        sliced_field_1 = getattr(sliced_p_1, f)
        sliced_field_2 = getattr(sliced_p_2, f)
        if field is not None:
            print(f"{p.shape}, {f}, {field.shape}")
            expected_sliced_shape_1 = field.shape[1:]
            expected_sliced_shape_2 = tuple([1] + list(field.shape[1:]))
            sliced_shape_1 = sliced_field_1.shape
            sliced_shape_2 = sliced_field_2.shape
            print(f"{f} one element: expected {expected_sliced_shape_1} vs got {sliced_shape_1}")
            print(f"{f} one slice: expected {expected_sliced_shape_2} vs got {sliced_shape_2}")
            assert sliced_shape_1 == expected_sliced_shape_1
            assert sliced_shape_2 == expected_sliced_shape_2
        else:
            print(f"{p.__class__.__name__}.{f} is None")
