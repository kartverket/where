"""Write the estimated state vector to screen.

Description:
------------

Write the final parameter estimates to screen. These estimates are corrections to the apriori values.
The printed number is just the mean value.




"""
# Standard library imports
from collections import OrderedDict

# External library imports
import numpy as np

# Where imports
from where.lib import log
from where.lib import plugins


@plugins.register
def parameter_corrections(dset):
    """Write statistics about baselines to file.

    Args:
        dset:   Dataset, information about model run.
    """
    state_vector_fields = sorted(
        [field for field in dset.fields if field.startswith("state_") and not field.endswith("_sigma")]
    )

    print_data = OrderedDict()
    for field in state_vector_fields:
        param_type, name = field[6:].split("-", maxsplit=1)
        name = name.split("_")
        field_sigma = field + "_sigma"
        key = (param_type, name[0], dset.unit(field))
        if len(name) == 1:
            print_data[key] = (np.mean(dset[field]), np.mean(dset[field_sigma]), 0)
        else:
            print_data.setdefault(key, dict())
            print_data[key][name[1]] = (
                np.mean(dset[field]),
                np.mean(dset[field_sigma]),
                np.sum(dset.filter(source=name[0])),
            )

    correction_str = ""
    for key, values in print_data.items():
        if isinstance(values, dict):
            suffix = sorted(values.keys())
            param = f"{key[0]}:{key[1]}_{''.join(suffix)}"
            correction_str += f"  {param:30} {values[suffix[0]][2]:4d} = "
            for k in suffix:
                correction_str += f"{values[k][0]:9.6f} "
            correction_str += f"{key[2]} (std:"
            for k in suffix:
                correction_str += f" {values[k][1]:6.4f}"
            correction_str += ")\n"
        else:
            correction_str += (
                f"  {f'{key[0]}:{key[1]}':30} {values[2]:4d} = {values[0]:9.6f} {key[2]} (std: {values[1]:6.4f})\n"
            )
    log.out(f"Parameter corrections from Kalman filter after estimation:\n{correction_str.rstrip()}")

    names = dset.meta["normal equation"]["names"]
    x = dset.meta["normal equation"]["solution"]
    Q = dset.meta["normal equation"]["covariance"]
    correction_str = ""
    for i, name in enumerate(names):
        correction_str += f"  {name:30} = {x[i]:15.12f} (std: {Q[i][i]:6.4f})\n"
    log.out(f"Parameter corrections from normal equations:\n{correction_str.rstrip()}")
