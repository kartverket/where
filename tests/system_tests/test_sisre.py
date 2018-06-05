"""Tests a full SISRE analysis

This test simply compares the result of a SISRE analysis with the result obtained earlier. In other words, it does NOT
guarantee that the results are correct. It will however test that the results have not changed.

"""

# System library imports
from datetime import date

# Third party imports
import pytest


@pytest.mark.sisre
@pytest.mark.system_test
@pytest.mark.slow
def test_sisre_analysis(system_test):
    """Test running a full SISRE analysis"""

    previous, current = system_test(
        rundate=date(2018, 2, 1), pipeline="sisre", stage="calculate", session="brdm", systems="E"
    )
    assert previous == current
