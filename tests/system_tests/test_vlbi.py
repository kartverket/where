"""Tests a full VLBI analysis

This test simply compares the result of a VLBI analysis with the result obtained earlier. In other words, it does NOT
guarantee that the results are correct. It will however test that the results have not changed.

"""

# System library imports
from datetime import date

# Third party imports
import pytest


@pytest.mark.vlbi
@pytest.mark.system_test
@pytest.mark.slow
def test_vlbi_analysis(system_test):
    """Test running a full VLBI analysis"""

    previous, current = system_test(rundate=date(2009, 11, 2), pipeline="vlbi", stage="estimate", session="XA")
    assert previous == current
