"""Basic functionality for parsing datafiles, extended by individual parsers

Description:

This module contains functions and classes for parsing datafiles in Where.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Standard library imports
from collections import UserDict
import itertools

# Where imports
from where.lib import config
from where.lib import dependencies
from where.lib import files
from where.lib import log
from where.lib.timer import timer
from where.lib import util


def define_parser(end_marker, label, parser_def, end_callback=None):
    """A convenience method for defining the necessary fields of a parser

    A single parser can read and parse one group of datalines, defined through this function by specifying how to parse
    each line (parser_def), how to identify each line (label), how to recognize the end of the group of lines
    (end_marker) and finally what (if anything) should be done after all lines in a group is read (end_callback).

    The end_marker, label and end_callback parameters should all be functions with the following signatures:

        end_marker   = func(line, line_num, next_line)
        label        = func(line, line_num)
        end_callback = func(cache)

    Args:
        end_marker:   A function returning True for the last line in a group.
        label:        A function returning a label used in the parser_def.
        parser_def:   A dict with 'parser' and 'fields' defining the parser.
        end_callback: A function called after reading all lines in a group.

    Returns:
        A dict containing the definition of the parser.
    """
    parser = dict(end_marker=end_marker, label=label, parser_def=parser_def)
    if end_callback is not None:
        parser["end_callback"] = end_callback

    return parser


class Parser(object):
    """An abstract base class that has basic methods for parsing a datafile

    This class provides functionality for parsing a file. You may also inherit from one of the Parser...-classes below.

    Attributes:
        file_key:       String, key to the datafile that will be read.
        rundate:        Datetime, the model run date.
        data_available: Indicator of whether data are available.
        vars:           Dict with variables used when identifying the datafile.
        data:           Dict containing the (observation) data read from file.
        meta:           Dict containing the metainformation read from file.
        dependencies:   List of files that have been read by the parser.

    """

    def __init__(self, rundate=None, file_path=None):
        """Set up the basic information needed by the parser

        Subclasses of Parser should extend this constructor, by calling super().__init__(rundate) before setting
        attributes (at the very least self.file_key).

        The `file_path` parameter is mainly for when running the parsers indepently of Where, and can be used to
        specify a file independent of `files.conf`. Note that this should not be done inside of the Where program, as
        that looses some of the logging and maintainability.

        Args:
            rundate (date):        The model run date (optional, used to set up date variables).
            file_path (String):    Optional path to file that will be read.
        """
        super().__init__()
        self.file_key = "Overwritten by subclasses"
        self.file_path = file_path
        self.rundate = rundate
        self.data_available = True
        self.dependencies = list()

        # Initialize the data
        self.vars = dict()
        self.meta = dict()
        self.data = dict()

        # Use _parser.Parser and subclasses instead
        log.dev(
            "parser.Parser is deprecated, let {} subclass one of LineParser, ChainParser or SinexParser instead",
            self.__class__.__name__,
        )

    def setup_parsers(self):
        util.not_implemented()

    def setup_calculators(self):
        return list()

    def parse(self):
        """Parse data

        This is a basic implementation that carries out the whole pipeline of reading and parsing datafiles including
        calculating secondary data.

        Returns:
            Parser: The parsed data
        """
        if self.file_path is None:
            self.file_path = files.path(self.file_key, file_vars=self.vars, download_missing=True)

        parser_package, parser_name = self.__module__.rsplit(".", maxsplit=1)
        with timer("Finish {} ({}) - {} in".format(parser_name, parser_package, self.file_key)):
            if self.data_available:
                self.read_data()

            if not self.data_available:  # May have been set to False by self.read_data()
                log.warn(
                    "No data found by {} for {} (was looking for {})",
                    self.__class__.__name__,
                    self.rundate.strftime(config.FMT_date),
                    self.file_path,
                )
                return self

            self.calculate_data()
            dependencies.add(*self.dependencies)

        return self

    def process_data(self):
        log.dev("{p}.process_data is deprecated. Use {p}.parse instead", p=self.__class__.__name__)
        self.parse()

    def read_data(self):
        """Read data from datafiles

        This is a basic implementation that parses data from one file as specified by the file_key and vars properties.

        For more advanced uses, e.g. reading data from several files, this method should be overridden. In that case,
        make sure file dependencies are appended to the self.dependencies-list.

        """
        self.dependencies.append(self.file_path)
        is_zipped = files.is_path_zipped(self.file_path)
        with files.open_path(self.file_path, mode="rt", is_zipped=is_zipped) as fid:
            self.parse_file(fid)

    def calculate_data(self):
        """
        TODO: Description?
        """
        for calculator in self.setup_calculators():
            log.debug("Start calculator {} in {}", calculator.__name__, self.__module__)
            with timer("Finish calculator {} ({}) in".format(calculator.__name__, self.__module__), logger=log.debug):
                calculator()

    def write_to_dataset(self, data_out):
        util.not_implemented()

    def parse_file(self, fid):
        """Read a data file and parse the contents
        """
        # Get "list" of parsers
        parsers_iter = iter(self.setup_parsers())

        parser = next(parsers_iter)  # Pointing to first parser
        cache = dict(line_num=0)

        # Get iterators for current and next line
        line_iter, next_line_iter = itertools.tee(fid)
        next(next_line_iter, None)

        # Iterate over all file lines including last line by using zip_longest
        for line, next_line in itertools.zip_longest(line_iter, next_line_iter):
            cache["line_num"] += 1
            self.parse_line(line.rstrip(), cache, parser)

            # Skip to next parser
            if next_line is None or parser["end_marker"](line.rstrip(), cache["line_num"], next_line):
                if "end_callback" in parser:
                    parser["end_callback"](cache)
                cache = dict(line_num=0)
                try:
                    parser = next(parsers_iter)
                except StopIteration:
                    break

    def parse_line(self, line, cache, parser):
        """Parse line

        A line is parsed by separating a line in fields. How the separation is done, is defined in the `parser_def`
        entry of the `parser` dictionary. The parser definition `parser_def` includes the `parser`, `field`, `strip` and
        `delimiter` entries. The `parser` entry points to the parser function and the `field` entry defines how to
        separate the line in fields. The separated fields are saved either in a dictionary or in a list. In the last
        case the line is splitted on whitespaces by default. With the `delimiter` entry the default definition can
        be overwritten. Leading and trailing whitespace characters are removed by default before a line is parsed.
        This default can be overwritten by defining the characters, which should be removed with the 'strip' entry. The
        `parser` dictionary is defined like:

          parser = { 'end_marker': <function>,
                     'label':      <function>,
                     'parser_def': { <label>: {'fields':    <dict or list of fields>,
                                               'parser':    <parser function>,
                                               'delimiter': <delimiter for splitting line>,
                                               'strip':     <characters to be removed from beginning or end of the line>
                                              }}}

        Args:
            line (str):     Line to be parsed.
            cache (dict):   Store temporary data.
            parser (dict):  Dictionary with defined parsers with the keys 'parser_def', 'label' and 'end_marker'.
        """
        if not parser["label"]:
            return

        # log.debug('{:>3d}: {}', cache['line_num'], line)

        label = parser["label"](line.rstrip(), cache["line_num"])
        if label not in parser["parser_def"]:
            return
        fields = parser["parser_def"][label]["fields"]
        values = dict()
        if isinstance(fields, dict):
            for field, idx in fields.items():
                values[field] = line[slice(*idx)].strip(parser["parser_def"][label].get("strip"))
        elif isinstance(fields, list):
            for field, value in zip(fields, line.split(parser["parser_def"][label].get("delimiter"))):
                if field is not None:
                    values[field] = value.strip(parser["parser_def"][label].get("strip"))

        parse_func = parser["parser_def"][label]["parser"]
        parse_func(values, cache)

    def parse_default(self, line, _):
        """Add the contents of line to data

        Args:
            line: Dict containing the fields of a line.
        """
        self.data.update(line)

    def parse_default_meta(self, line, _):
        """Add the contents of line to meta

        Args:
            line: Dict containing the fields of a line.
        """
        self.meta.update(line)

    def copy_cache_to_data(self, cache):
        """Copy contents of the cache to the data datastructure

        The fields 'key' and 'values' must exist in the cache. The value of cache['key'] is used as a field name in
        data, to which the contents of cache['values'] are added. The value of cache['values'] should be a dictionary.

        Args:
            cache:   Dictionary with the fields 'key' and 'values'.

        """
        self.data.setdefault(cache["key"], {}).update(cache["values"])

    def copy_cache_to_meta(self, cache):
        """Copy contents of the cache to the meta datastructure

        The fields 'key' and 'values' must exist in the cache. The value of cache['key'] is used as a field name in
        meta, to which the contents of cache['values'] are added. The value of cache['values'] should be a dictionary.

        Args:
            cache:   Dictionary with the fields 'key' and 'values'.

        """
        self.meta.setdefault(cache["key"], {}).update(cache["values"])

    def __hash__(self):
        return id(self)


class ParserDict(Parser, UserDict):
    """A Parser that can be treated as a dict containing all the parsed data
    """
    pass
