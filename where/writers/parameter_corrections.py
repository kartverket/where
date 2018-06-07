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
                np.mean(dset[field]), np.mean(dset[field_sigma]), np.sum(dset.filter(source=name[0]))
            )

    print("Parameter corrections from Kalman filter after estimation:")
    for key, values in print_data.items():
        if isinstance(values, dict):
            suffix = sorted(values.keys())
            print("{:30} {:4d} = ".format(key[0] + ":" + key[1] + "_" + "".join(suffix), values[suffix[0]][2]), end="")
            for k in suffix:
                print("{: 14.12f} ".format(values[k][0]), end="")
            print(key[2], "(std:", end="")
            for k in suffix:
                print("{: 6.4f} ".format(values[k][1]), end="")
            print(")")
        else:
            print(
                "{:30} {:4d} = {: 14.12f} {} (std: {:6.4f})".format(
                    key[0] + ":" + key[1], values[2], values[0], key[2], values[1]
                )
            )

    print("Parameter corrections from normal equations")
    names = dset.meta["normal equation"]["names"]
    x = dset.meta["normal equation"]["solution"]
    Q = dset.meta["normal equation"]["covariance"]

    for i, name in enumerate(names):
        print("{:30} = {: 14.12f} (std: {: 6.4f})".format(name, x[i], Q[i][i]))
