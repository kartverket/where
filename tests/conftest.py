"""Common functions for all tests

"""

# System library imports
from datetime import datetime
import shutil

# Third party imports
import pytest

# Where imports
from where.lib import config


@pytest.fixture
def backup_test_file():
    """Copy a given file to a work/test-subdirectory"""
    return _backup_test_file


def _backup_test_file(source_path, subdirectory):
    """Copy a given file to a work/test-subdirectory

    Args:
        source_path (Path):   Path to file that will be copied.
        subdirectory (str):   Subdirectory beneath {path_work}/test where the file will be copied.
    """
    print("BTF", source_path, subdirectory)

    # Create a unique name with a timestamp
    destination_dir = config.where.path.work.path / "test" / subdirectory
    destination_name = source_path.name + f".{datetime.now():%Y%m%d-%H%M%S}"
    destination_path = destination_dir / destination_name

    # Copy source to destination
    destination_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(source_path, destination_path)
