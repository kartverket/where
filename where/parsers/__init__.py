"""Framework for parsers

Description:
------------

To add a new parser, simply create a new .py-file which defines a class inheriting from parsers.Parser. The class needs
to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.parsers import parser
    from where.lib import plugins

    @plugins.register
    class MyNewParser(parser.Parser):
        ...

To use a parser, you will typically use the :func:`parse`-function defined below

    from where import parsers
    my_new_parser = parsers.parse('my_new_parser', ...)
    my_data = my_new_parser.as_dict()

The name used in `parse` to call the parser is the name of the module (file) containing the parser.




$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# Where imports
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib import plugins


def names():
    """List the names of the available parsers

    Returns:
        List: List of strings with the names of the available parsers
    """
    return plugins.list_all(package_name=__name__)


def setup_parser(parser_name=None, file_key=None, **kwargs):
    """Set up the given parser.

    Note that this only sets up the parser, no data will be read and parsed.

    It is possible to give file key instead of parser name. In that case the name of the parser will be read from the
    file list.

    Args:
        parser_name (String):   Name of parser.
        file_key (String):      Used to look up parser in the Where file list.
        kwargs:                 Input arguments to the parser.

    Returns:
        Parser:  An instance of the given parser

    """
    parser_name = config.files.get(section=file_key, key="parser", value=parser_name).str
    log.assert_not_none(parser_name, "No parser found for '{}' in {}", file_key, ", ".join(config.files.sources))

    parser = plugins.call_one(package_name=__name__, plugin_name=parser_name, use_timer=False, **kwargs)

    if file_key is not None:
        parser.file_key = file_key

    return parser


def parse(parser_name=None, file_key=None, **kwargs):
    """Call the given parser and return parsed data

    It is possible to give file key instead of parser name. In that case the name of the parser will be read from the
    file list.

    Args:
        parser_name (String):   Name of parser
        file_key (String):      Used to look up parser in the Where file list.
        kwargs:                 Input arguments to the parser

    Returns:
        Parser:  The parsed data
    """
    return setup_parser(parser_name=parser_name, file_key=file_key, **kwargs).parse()


def parse_file(parser_name, file_path, use_cache=True, **parser_args):
    """Use the given parser on a file and return parsed data

    Specify `parser_name` and `file_path` to the file that should be parsed. The following parsers are available:

    {doc_parser_names}

    Data can be retrieved either as Dictionaries, Pandas DataFrames or Where Datasets by using one of the methods
    `as_dict`, `as_dataframe` or `as_dataset`.

    Example:
        > df = parsers.parse_file('rinex_obs', 'ande3160.16o').as_dataframe()

    Args:
        parser_name (String):     Name of parser
        file_path (String/Path):  Path to file that should be parsed.
        use_cache (Boolean):      Whether to use a cache to avoid parsing the same file several times.
        parser_args:              Input arguments to the parser

    Returns:
        Parser:  Parser with the parsed data
    """
    # Create the parser and parse the data
    parser = plugins.call_one(
        package_name=__name__, plugin_name=parser_name, use_timer=False, file_path=file_path, **parser_args
    )
    return parser.parse()


def parse_key(file_key, file_vars=None, parser=None, use_cache=True, **parser_args):
    """Parse a file given in the Where file-list and return parsed data

    By specifying a `file_key`. The file_key is looked up in the file list to figure out which file that should be
    parsed. The name of the parser will also be looked up in the file configuration. The dictionary `file_vars` may be
    specified if variables are needed to figure out the correct file path from the configuration. The following file
    keys are available:

    {doc_file_keys}

    Data can be retrieved either as Dictionaries, Pandas DataFrames or Where Datasets by using one of the methods
    `as_dict`, `as_dataframe` or `as_dataset`.

    Example:
        > df = parsers.parse_key('center_of_mass', file_vars=dict(satellite='Lageos')).as_dataset()

    Args:
        file_key (String):        Used to look up parser_name and file_path in the Where file configuration.
        file_vars (Dict):         Additional file variables used when looking up file path in configuration.
        use_cache (Boolean):      Whether to use a cache to avoid parsing the same file several times.
        parser_args:              Input arguments to the parser

    Returns:
        Parser:  Parser with the parsed data
    """
    # Read parser_name from config.files if it is not given
    parser_name = config.files.get(section=file_key, key="parser", value=parser).str
    log.assert_not_none(parser_name, "No parser found for '{}' in {}", file_key, ", ".join(config.files.sources))

    # Figure out the file path
    file_vars = dict() if file_vars is None else file_vars
    file_path = files.path(file_key, file_vars=file_vars, download_missing=True, use_aliases=True)

    # Create the parser and parse the data
    return parse_file(parser_name, file_path, use_cache=use_cache, **parser_args)
