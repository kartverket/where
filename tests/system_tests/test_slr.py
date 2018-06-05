"""Tests a full SLR analysis

This test simply compares the result of a SLR analysis with the result obtained earlier. In other words, it does NOT
guarantee that the results are correct. It will however test that the results have not changed.

"""

# System library imports
from datetime import date

# Third party imports
import pytest


@pytest.mark.slr
@pytest.mark.system_test
@pytest.mark.slow
def test_slr_analysis(system_test):
    """Test running a full SLR analysis"""

    previous, current = system_test(
        rundate=date(2015, 9, 1), pipeline="slr", stage="calculate", session="lageos1", arc_length=1
    )
    assert previous == current
