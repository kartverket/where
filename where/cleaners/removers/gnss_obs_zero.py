"""Remove large observation outliers

Description:

Todo

"""
# Standard library imports

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import plugins


@plugins.register
def gnss_obs_zero(dset):
    """Edits data based on data quality

    NOTE: This only a workaround and works only if observations from one GNSS are used. If several GNSSs are used then
          it can happen that due to the fact that some observation type are only given for a certain GNSS, that all
          observations are deleted.

    Args:
        dset:     A Dataset containing model data.

    Returns:
        Numpy array:  False for observations to throw away.
    """
    session = dset.dataset_name
    flag = config.session[session].gnss_obs_zero.bool
    edit_dset = np.full(dset.num_obs, True, dtype=bool)

    if flag is False:  # TODO: Should it be done like that, if editor is not in use?
        return edit_dset

    # Loop over GNSSs and observation types
    for sys in dset.meta["obstypes"]:
        for obstype in dset.meta["obstypes"][sys]:

            # Exclude Doppler and SNR observations
            if (obstype not in dset.fields) or (obstype[0] in "DS"):
                continue

            # Remove observations with (close to) zero value
            edit_dset = np.logical_and(edit_dset, dset[obstype] > 10000)

    return edit_dset
