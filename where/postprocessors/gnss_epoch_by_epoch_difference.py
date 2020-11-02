"""Add epoch by epoch difference of observations to dataset

Description:
------------

"""
# External imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_epoch_by_epoch_difference(dset: "Dataset") -> None:
    """Add epoch by epoch difference of observations to dataset

    Args:
        dset:     A Dataset containing model data.
    """
    for field in dset.obs.fields:

        log.debug(f"Add epoch by epoch difference for observation field '{field}' to dataset.")        
        diff = np.full(dset.num_obs, float('nan'))
        
        for sys in dset.unique("system"):      
            idx_sys = dset.filter(system=sys)
            sys_tmp = dset.obs[field][idx_sys]
                      
            for sat in set(dset.satellite[idx_sys]):          
                idx_sat = dset.satellite[idx_sys] == sat
                sys_tmp[idx_sat] = np.insert(
                                        np.diff(dset.obs[field][idx_sys][idx_sat]), 
                                        0, 
                                        float('nan'),
                ) # first epoch is not differenced, therefore NaN has to be inserted as first element.
                
            diff[idx_sys] = sys_tmp
            
        dset.add_float(f"diff_epo.{field}", val=diff, unit=dset.unit(f"obs.{field}"))
        

            