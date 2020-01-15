"""Common functions for system tests

"""

# System library imports
import subprocess

# Third party imports
import pytest

# Where imports
import where
from where.lib import config


@pytest.fixture
def where_runner():
    return _where_runner


def _where_runner(command):
    """Run a session of Where
    """
    print("RUN", command)
    subprocess.run(command.split())


@pytest.fixture
def parameters():
    """Set up parameters for running Where"""
    return _parameters


def _parameters(rundate, pipeline, stage, session, **options):
    """Set up parameters for running Where

    Args:
        rundate (date):   The model run date.
        pipeline (str):   The pipeline.
        stage (str):      The stage to compare.
        session (str):    The session to compare.
        options (dict):   Command line options that will be passed to Where.

    Returns:
        dict: Command and file_vars representing the given Where analysis
    """
    user = "test"
    params = dict(
        command=(
            f"{where.__executable__} {rundate:%Y %m %d} --{pipeline} --session={session} -D -N "
            f"--user={user} --output=system_test:{stage} " + " ".join(f"--{o}={v}" for o, v in options.items())
        ),
        file_vars=dict(
            rundate=rundate.strftime(config.FMT_date),
            pipeline=pipeline,
            stage=stage,
            user=user,
            id="",
            session=session,
            **config.date_vars(rundate),
        ),
    )

    return params


@pytest.fixture
def system_test(backup_test_file):
    """Perform a system test by running a complete analysis in Where, and comparing the result with last analysis"""

    def _system_test(rundate, pipeline, stage, session, **options):
        """Perform a system test

        Args:
            rundate (date):   The model run date.
            pipeline (str):   The pipeline.
            stage (str):      The stage to compare.
            session (str):    The session to compare.
            options (dict):   Command line options that will be passed to Where.

        Returns:
            tuple: Previous and Current output results.
        """
        parameters = _parameters(rundate, pipeline, stage, session, **options)

        # File containing output result
        file_path = config.files.path("output_system_test", file_vars=parameters["file_vars"])

        # Read previous result
        try:
            previous_result = file_path.read_text()
            if not previous_result:
                previous_result = "Empty file"
        except FileNotFoundError:
            previous_result = "No prior result found"

        # Run analysis
        _where_runner(parameters["command"])

        # Keep a backup of the result file
        backup_test_file(file_path, "system_test")

        # Compare result
        current_result = file_path.read_text()

        # Workaround for pytest sometimes hanging when comparing long strings
        for prev, cur in zip(previous_result.split("\n"), current_result.split("\n")):
            if prev != cur:
                return prev, cur

        return previous_result, current_result

    return _system_test
