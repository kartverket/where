"""Basic functionality for parsing datafiles line by line using Numpy

Description:

This module contains functions and classes for parsing datafiles in Where.


"""
# Standard library imports
from datetime import datetime

# External library imports
import numpy as np

# Where imports
from where.parsers._parser import Parser
from where.lib import util


class LineParser(Parser):
    """An abstract base class that has basic methods for parsing a datafile

    This class provides functionality for using numpy to parse a file line by line. You should inherit from this one,
    and at least specify the necessary parameters in `setup_parser`.
    """

    def __init__(self, file_path, encoding=None):
        """Set up the basic information needed by the parser

        Add a self._array property for the raw numpy array data.

        Args:
            file_path (String/Path):    Path to file that will be read.
        """
        super().__init__(file_path, encoding)
        self._array = None

    def setup_parser(self):
        """Set up information needed for the parser

        This should return a dictionary with all parameters needed by np.genfromtxt to do the actual parsing.

        Returns:
            Dict:  Parameters needed by np.genfromtxt to parse the input file.
        """
        util.not_implemented()

    def read_data(self):
        """Read data from the data file

        Uses the np.genfromtxt-function to parse the file. Any necessary parameters should be returned by
        `setup_parser`. See `self.structure_data` if the self.data-dictionary needs to be structured in a particular
        way.
        """
        self.meta["__params__"] = self.setup_parser()
        self.meta["__params__"].setdefault("encoding", self.file_encoding or "bytes")  # TODO: Default to None instead?
        self._array = np.genfromtxt(self.file_path, **self.meta["__params__"])
        self.structure_data()

    def structure_data(self):
        """Structure raw array data into the self.data dictionary

        This simple implementation creates a dictionary with one item per column in the array. Override this method for
        more complex use cases.
        """
        for name in self._array.dtype.names:
            self.data[name] = self._array[name]
