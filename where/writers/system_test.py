"""Write output that can be used in system tests

Description:
------------

Write simple output from an analysis.


"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import files


@plugins.register
def system_test_output(dset):
    """Write simple output based on dataset

    Args:
        dset:   Dataset, information about model run.
    """
    fields = config.tech.get("fields", section="system_test").tuple

    with files.open("output_system_test", file_vars=dset.vars, mode="wt") as fid:
        for idx, vals in enumerate(dset.values(*fields), start=1):
            fid.write(f"{idx:6d} " + " ".join(str(v) for v in vals) + "\n")
