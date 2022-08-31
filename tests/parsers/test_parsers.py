"""Tests for the parsers-package

Tests for individual parsers are in their own files

Example:
--------
    python -m pytest -s test_parsers.py
"""

# Standard library imports
from datetime import datetime
import pathlib

# Where imports
from where import parsers
from where.lib import config


#
# TEST CONFIGURATION
#

# Initialize configuration
#
# TODO: Configuration handling should be improved. Use a test_parsers.conf file.
config.init(
    rundate = datetime(2021, 2, 1),
    pipeline = "gnss",
)


def get_parser(parser_name, example_path = None):
    """Get a parser that has parsed an example file"""
    if not example_path:
        example_path = pathlib.Path(__file__).parent / "example_files" / parser_name
    return parsers.parse_file(parser_name, example_path)


#
# Tests
#
def test_list_of_parsers():
    """Test that names of parsers can be listed"""
    assert len(parsers.names()) > 0



def test_parser_gnss_has_decoder():
    """Test that parsing gnss_has_decoder gives expected output"""
    parser = get_parser(
        "gnss_has_decoder", 
        example_path=pathlib.Path(__file__).parent / "example_files" / "SEPT147.21_has_cb.csv",
    ).as_dict()

    assert len(parser) == 13
    assert "ToW" in parser
    assert 259 in parser["ToW"]
