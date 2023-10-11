"""Apply Galileo HAS corrections to RINEX observations

"""

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import config, log, util

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def rinex_obs_apply_has_correction(dset: "Dataset") -> None:
    """Apply Galileo HAS corrections to RINEX observations

    Args:
        dset (Dataset):       A dataset containing the data.
    """
    
    # Apply HAS code and phase bias correction        
    file_path = config.files.path("gnss_has_cb", file_vars={**dset.vars, **dset.analysis})
    if file_path.exists() and not util.is_file_empty(file_path):
        code_bias_has = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            file_key="gnss_has_cb",
            day_offset=0, #MURKS: Should be 1
            apriori_orbit="has",
        )
        code_bias_has.apply_code_bias_to_dataset(dset)
    else:
        log.warn(f"File path {file_path} does not exists or is empty. That means no code bias HAS corrections are "
                 "applied.")

    file_path = config.files.path("gnss_has_cp", file_vars={**dset.vars, **dset.analysis})
    if file_path.exists() and not util.is_file_empty(file_path):
        phase_bias_has = apriori.get(
            "orbit", 
            rundate=dset.analysis["rundate"], 
            file_key="gnss_has_cp",
            day_offset=0, #MURKS: Should be 1
            apriori_orbit="has",
        )
        phase_bias_has.apply_phase_bias_to_dataset(dset)
    else:
        log.warn(f"File path {file_path} does not exists or is empty. That means no phase bias HAS corrections are "
                 "applied.")

