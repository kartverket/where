"""Definition of Where-specific exceptions

Description:
------------

Custom exceptions used by Where for more specific error messages and handling.

"""

from midgard.dev.exceptions import MidgardException  # noqa


class WhereException(Exception):
    pass


class WhereExit(SystemExit, WhereException):
    pass


class InitializationError(WhereException):
    pass


class InvalidArgsError(WhereException):
    pass


class FieldExistsError(WhereException):
    pass


class FieldDoesNotExistError(WhereException):
    pass


class MissingDataError(WhereException):
    pass


class UnknownEnumError(WhereException):
    pass


class UnknownPipelineError(WhereException):
    pass


class UnknownRadioSourceError(WhereException):
    pass


class UnknownSiteError(WhereException):
    pass


class UnitError(WhereException):
    pass
