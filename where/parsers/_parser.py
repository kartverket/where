"""Basic functionality for parsing datafiles, extended by individual parsers

Description:

This module contains functions and classes for parsing datafiles in Where.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""
# Standard library imports
import numbers
import pathlib

# External library imports
import pandas as pd

# Where imports
from where.lib import dependencies
from where.lib import log
from where.lib.timer import timer
from where.lib import util


class Parser:
    """An abstract base class that has basic methods for parsing a datafile

    This class provides functionality for parsing a file. You should inherit from one of the specific parsers like for
    instance ChainParser, LineParser, SinexParser etc

    Attributes:
        file_path (Path/String):      Path to the datafile that will be read.
        parser_name (String):         Name of the parser (as needed to call parsers.parse_...).
        data_available (Boolean):     Indicator of whether data are available.
        data (Dict):                  The (observation) data read from file.
        meta (Dict):                  Metainformation read from file.
    """

    def __init__(self, file_path):
        """Set up the basic information needed by the parser

        Args:
            file_path (String/Path):    Path to file that will be read.
        """
        self.file_path = pathlib.Path(file_path)
        self.parser_name = self.__module__.split(".")[-1]

        # Initialize the data
        self.data_available = self.file_path.exists()
        self.meta = dict(__parser_name__=self.parser_name, __data_path__=self.file_path)
        self.data = dict()

    def setup_parser(self):
        pass

    def setup_calculators(self):
        return list()

    def parse(self):
        """Parse data

        This is a basic implementation that carries out the whole pipeline of reading and parsing datafiles including
        calculating secondary data.

        Subclasses should typically implement (at least) the `read_data`-method.
        """
        self.setup_parser()

        parser_package = self.__module__.rsplit(".", maxsplit=1)[0]
        with timer("Finish {} ({}) - {} in".format(self.parser_name, parser_package, self.file_path)):
            if self.data_available:
                self.read_data()

            if not self.data_available:  # May have been set to False by self.read_data()
                log.warn("No data found by {} in {}", self.__class__.__name__, self.file_path)
                return self

            self.calculate_data()

        dependencies.add(self.file_path)
        return self

    def read_data(self):
        """Read data from the data file

        Data should be read from `self.file_path` and stored in the dictionary `self.data`. A description of the data
        may be placed in the dictionary `self.meta`. If data are not available for some reason, `self.data_available`
        shoulde be set to False.
        """
        util.not_implemented()

    def calculate_data(self):
        """Do simple manipulations on the data after they are read

        Simple manipulations of data may be performed in calculators after they are read. They should be kept simple so
        that a parser returns as true representation of the data file as possible. Advanced calculations may be done
        inside apriori classes or similar.

        To add a calculator, define it in its own method, and override the `setup_calculators`-method to return a list
        of all calculators.
        """
        for calculator in self.setup_calculators():
            log.debug("Start calculator {} in {}", calculator.__name__, self.__module__)
            with timer("Finish calculator {} ({}) in".format(calculator.__name__, self.__module__), logger=log.debug):
                calculator()

    def as_dict(self, include_meta=False):
        """Return the parsed data as a dictionary

        Args:
            include_meta (Boolean):   Whether to include meta-data in the returned dictionary (default: False).

        Returns:
            Dictionary:  The parsed data.
        """
        return dict(self.data, __meta__=self.meta) if include_meta else self.data.copy()

    def as_dataframe(self, index=None):
        """Return the parsed data as a Pandas DataFrame

        This is a basic implementation, assuming the `self.data`-dictionary has a simple structure. More advanced
        parsers may need to reimplement this method.

        Args:
            index (String / List):      Name of field to use as index. May also be a list of strings.

        Returns:
            DataFrame: The parsed data.
        """
        df = pd.DataFrame.from_dict(self.data)
        if index is not None:
            df.set_index(index, drop=True, inplace=True)

        return df

    def as_dataset(self):
        """Return the parsed data as a Where Dataset

        This is a basic implementation, assuming the `self.data`-dictionary has a simple structure. More advanced
        parsers may need to reimplement this method.

        Returns:
            Dataset:  The parsed data.
        """
        from where import data

        num_obs = len(self.data[list(self.data.keys()).pop()])  # Take num obs from one random field
        dset = data.Dataset.anonymous(num_obs=num_obs)
        for field in self.data.keys():
            if isinstance(self.data[field][0], numbers.Number):
                dset.add_float(field, val=self.data[field])
            elif isinstance(self.data[field][0], str):
                dset.add_text(field, val=self.data[field])
        return dset

    def update_dataset(self, dset):
        """Update the given dataset with the parsed data

        This is a basic implementation, assuming the `self.data`-dictionary has a simple structure. More advanced
        parsers may need to reimplement this method.

        Args:
            dset (Dataset):  The dataset to update with parsed data.
        """
        # parser_dset = self.as_dataset()
        # if new fields:
        #     dset.add ...
        # elif new epochs:
        #     dset.extend ...
        util.not_implemented()

    def __repr__(self):
        return "{}(file_path='{}')".format(self.__class__.__name__, self.file_path)
