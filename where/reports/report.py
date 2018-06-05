"""WHERE library module for storing data for a report

Description:
------------

asdf



$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""

# Standard library imports
import atexit
from collections import OrderedDict
from datetime import datetime
import pickle

# WHERE imports
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib import util

# _REPORTS is a list of reports
_REPORTS = OrderedDict()


def add(report_name, dset=None, **report_vars):
    """Add general report to report list

    The message is parsed with str.format and written to the report.

    Args:
        report_type:   Message log level, one of the keys of #LOGLEVELS.
        report_text:   String, will be parsed with str.format.
        formatargs:    Arguments passed on to str.format.
        formatkws:     Keyword arguments passed on to str.format.
    """
    # Do not add anything if report.init has not been called
    if "__section__" not in _REPORTS:
        return

    # Add report to current section
    section = _REPORTS["__section__"]
    log.debug("Adding {} to {} for report", report_name, section)
    _REPORTS[section].append(
        dict(report_name=report_name, dset_params={} if dset is None else dset.parameters, **report_vars)
    )


def data(report_name="data", dset=None, **data):
    add(report_name, dset=dset, **data)


def text(report_name="text", text="", *format_args, **format_kwargs):
    formatted_text = text.format(*format_args, **format_kwargs)
    add(report_name, text=formatted_text)


def start_session(session):
    _REPORTS["__session__"] = session
    start_section("DEFAULT")


def start_section(section):
    full_section_name = "{session}:{section}".format(session=_REPORTS["__session__"], section=section)
    _REPORTS["__section__"] = full_section_name
    _REPORTS[full_section_name] = list()


def init(sessions):
    """Start the reporting.
    """
    if not config.tech.report.str:
        return

    read_reports_from_file()
    atexit.register(write_reports_to_file())

    rundate = config.analysis.rundate.date
    tech = config.analysis.tech.str

    # Remove sessions that will not be used
    keys_to_delete = [k for k in _REPORTS.keys() if ":" in k and k.split(":")[0] not in sessions]
    for key in keys_to_delete:
        del _REPORTS[key]

    # Write reports to DEFAULT until session or section is specified
    start_session("DEFAULT")

    # Should this be more general?
    text("title", "Where {:%Y-%m-%d} for {}", rundate, tech.upper())
    text(
        "text",
        "The model run was started at {} using the command `{}`.",
        datetime.now().strftime(config.FMT_datetime),
        util.COMMAND,
    )
    text(
        "code",
        "# Basic set up\n%matplotlib inline\n"
        "import where\nfrom where import models\nfrom where.lib import dataset\n"
        "from datetime import datetime\n"
        "import numpy as np\nimport matplotlib.pyplot as plt",
    )

    text("header", "Configuration")
    text("code", "import where\nwhere.set_config({rd:%Y}, {rd:%m}, {rd:%d}, '{tech}')", rd=rundate, tech=tech)

    text("header", "Run full model")
    text("text", "The following code will run the full model. Note that this might overwrite this notebook.")
    text("code", "# where.run({rd:%Y}, {rd:%m}, {rd:%d}, '{tech}', '-F')", rd=rundate, tech=tech)


def read_reports_from_file():
    try:
        with files.open("report_pickle", mode="rb") as fid:
            tmp_reports = pickle.load(fid)
            for section, values in tmp_reports.items():
                if section == "DEFAULT" or section.startswith("_"):
                    continue
                _REPORTS[section] = values

    except FileNotFoundError:
        pass


def write_reports_to_file():
    # Store current config in case it changes before write_to_file-function is actually run
    tech = config.analysis.tech.str
    rundate = config.analysis.rundate.str
    file_vars = config.files.vars.copy()

    def write_to_file():
        log.debug("Store reports for {} {}".format(tech.upper(), rundate))
        with files.open("report_pickle", mode="wb", file_vars=file_vars) as fid:
            pickle.dump(_REPORTS, fid)

    return write_to_file


def reports(*filter_types):
    from where.data import Dataset

    for section in _REPORTS:
        if section.startswith("_"):
            continue

        for report in _REPORTS[section]:
            name = "{}/{}".format(section, report["report_name"])
            if not filter_types or report["report_name"] in filter_types or name in filter_types:
                report_dset = Dataset(**report["dset_params"]) if report["dset_params"] else None
                yield name, report_dset, report
