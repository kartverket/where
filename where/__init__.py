"""Where, analysis of space geodetic data

This module provides for interactive use of the Where program. For running the Where program as a command line script,
see :mod:`where.__main__`.

Current Maintainers:
--------------------

{maintainers}

"""

# Standard library imports
from datetime import date as _date
from collections import namedtuple as _namedtuple
import platform as _platform


# Version of Where.
#
# This is automatically set using the bumpversion tool
__version__ = "2.1.2"


# Authors of the software
_Author = _namedtuple("_Author", ["name", "email", "start", "end"])

_AUTHORS = [
    _Author("Michael Dähnn", "michael.daehnn@kartverket.no", _date.min, _date.max),
    _Author("Ingrid Fausk", "ingrid.fausk@kartverket.no", _date.min, _date.max),
    _Author("Ann-Silje Kirkvik", "ann-silje.kirkvik@kartverket.no", _date.min, _date.max),
    _Author("Mohammed Ouassou", "mohammed.ouassou@kartverket.no", _date(2018, 9, 1), _date.max),
    # Hall of Fame
    _Author("Eirik Mysen", "eirik.mysen@kartverket.no", _date.min, _date(2017, 6, 1)),
    _Author("Geir Arne Hjelle", "geir.arne.hjelle@kartverket.no", _date.min, _date(2019, 2, 1)),
]

__author__ = ", ".join(a.name for a in _AUTHORS if a.start < _date.today() < a.end)
__contact__ = ", ".join(a.email for a in _AUTHORS if a.start < _date.today() < a.end)


# Copyleft of the software
__copyright__ = "2015 - {} Kartverket".format(_date.today().year)


# Name of executable (on Windows/DOS, WHERE is a builtin command), add custom formatter to list subprograms
class Executable(str):

    _executables = {
        "where": "{exe}",
        "there": "{gui_exe}",
        "profiler": "{exe}_profiler",
        "release": "{exe}_release",
        "runner": "{exe}_runner",
        "tools": "{exe}_tools",
    }

    def __new__(cls, executable, gui_executable):
        self = super().__new__(cls, executable)
        self.executables = {"exe": executable, "gui_exe": gui_executable}
        return self

    def __format__(self, format_spec):
        if not format_spec:
            format_spec = "where"

        return self._executables.get(format_spec, "").format(**self.executables)


if _platform.system() == "Windows":
    __executable__ = Executable("gd_where", "gd_there")
else:
    __executable__ = Executable("where", "there")


# Update doc with info about maintainers
def _update_doc(doc):
    """Add information to doc-string

    Args:
        doc (str):  The doc-string to update.

    Returns:
        str: The updated doc-string
    """
    # Maintainers
    maintainer_list = [f"+ {a.name} <{a.email}>" for a in _AUTHORS if a.start < _date.today() < a.end]
    maintainers = "\n".join(maintainer_list)

    # Add to doc-string
    return doc.format(maintainers=maintainers)


__doc__ = _update_doc(__doc__)
