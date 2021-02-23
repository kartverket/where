"""Framework for parsers

Description:
------------

To add a new parser, simply create a new .py-file which defines a class inheriting from parsers.Parser. The class needs
to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from where.parsers import _parser
    from midgard.dev import plugins

    @plugins.register
    class MyNewParser(parser.Parser):
        ...

To use a parser, you will typically use the :func:`parse_key`- or :func:`parse_file-functions defined below

    from where import parsers
    my_new_parser = parsers.parse_file('my_new_parser', 'file_name.txt', ...)
    my_data = my_new_parser.as_dict()

The name used in `parse_file` to call the parser is the name of the module (file) containing the parser.

"""

# Midgard imports
from midgard.parsers import names, parse_file  # noqa
from midgard import parsers as mg_parsers
from midgard.dev import plugins
from midgard.files import dependencies

# Where imports
from where.lib import config
from where.lib import log

# Add Where parsers to Midgard parsers
plugins.add_alias(mg_parsers.__name__, __name__)


def setup_parser(parser_name=None, file_key=None, **kwargs):
    """Set up the given parser.

    Note that this only sets up the parser, no data will be read and parsed.

    It is possible to give file key instead of parser name. In that case the name of the parser will be read from the
    file list.

    TODO: This is the old style of running parsers, can be deleted when all parsers are new style.

    Args:
        parser_name (String):   Name of parser.
        file_key (String):      Used to look up parser in the Where file list.
        kwargs:                 Input arguments to the parser.

    Returns:
        Parser:  An instance of the given parser

    """
    parser_name = config.files.get(section=file_key, key="parser", value=parser_name).str
    if not parser_name:
        log.warn(f"No parser found for {file_key!r} in {', '.join(config.files.sources)}")

    parser = plugins.call(package_name=mg_parsers.__name__, plugin_name=parser_name, **kwargs)

    if file_key is not None:
        parser.file_key = file_key

    return parser


def parse(parser_name=None, file_key=None, **kwargs):
    """Call the given parser and return parsed data

    It is possible to give file key instead of parser name. In that case the name of the parser will be read from the
    file list.

    TODO: This is the old style of running parsers, can be deleted when all parsers are new style.

    Args:
        parser_name (String):   Name of parser
        file_key (String):      Used to look up parser in the Where file list.
        kwargs:                 Input arguments to the parser

    Returns:
        Parser:  The parsed data
    """
    return setup_parser(parser_name=parser_name, file_key=file_key, **kwargs).parse()


def parse_key(file_key, file_vars=None, parser_name=None, use_cache=True, **parser_args):
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
        file_key (String):     Used to look up parser_name and file_path in the Where file configuration.
        file_vars (Dict):      Additional file variables used when looking up file path in configuration.
        parser_name (String):  Name of parser to use. Default is to use parser named in the file list.
        use_cache (Boolean):   Whether to use a cache to avoid parsing the same file several times.
        parser_args:           Input arguments to the parser.

    Returns:
        Parser:  Parser with the parsed data
    """
    # Read parser_name from config.files if it is not given
    parser_name = config.files.get(section=file_key, key="parser", value=parser_name).str
    if not parser_name:
        log.warn(f"No parser found for {file_key!r} in {', '.join(config.files.sources)}")

    # Figure out the file path
    file_vars = dict() if file_vars is None else file_vars
    download_missing = config.where.files.download_missing.bool
    file_path = config.files.path(file_key, file_vars=file_vars, download_missing=download_missing, use_aliases=True)
    dependencies.add(file_path, label=file_key)
    parser_args.setdefault("encoding", config.files.encoding(file_key))

    # Use the Midgard parser function to create parser and parse data
    return parse_file(parser_name, file_path, use_cache=use_cache, timer_logger=log.time, **parser_args)
