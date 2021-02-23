"""Add GNSS linar observation combinations to dataset

Description:
------------

Depending on the configuration, following linear combination can be added:
    code_multipath
    code_phase
    geometry_free
    ionosphere_free
    melbourne_wuebbena
    narrow_lane
    wide_lane
    
"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.lib.gnss import code_phase_difference, linear_combination, linear_combination_cmc, linear_combination_melbourne

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def gnss_linear_combination(dset: "Dataset") -> None:
    """Add GNSS linar observation combinations to dataset

    Args:
        dset:     A Dataset containing model data.
    """
    
    func = {
        "code_multipath": linear_combination_cmc,
        "code_phase": code_phase_difference,
        "geometry_free": linear_combination,
        "ionosphere_free": linear_combination,
        "melbourne_wuebbena": linear_combination_melbourne,
        "narrow_lane": linear_combination,
        "wide_lane": linear_combination,
    }
    
    for comb_name in config.tech[_SECTION].linear_combination.list:
        
        log.debug(f"Add {comb_name} combination to dataset.")
        
        # Code-multipath linear combination
        if comb_name == "code_multipath":
            try:
                cmc1, cmc2 = func[comb_name](dset)
            except ValueError:
                log.warn(f"Code multipath linear combination is not added to dataset. Dual-frequency code and phase "
                         f"observations are needed.")
                continue 
            
            dset.add_float(f"lin.{comb_name}_f1", val=cmc1["val"], unit="meter")
            dset.add_float(f"lin.{comb_name}_f2", val=cmc2["val"], unit="meter")
            dset.meta.setdefault("linear_combination", dict()).update({f"{comb_name}_f1": cmc1["sys_obs"]})
            dset.meta["linear_combination"][f"{comb_name}_f2"] = cmc2["sys_obs"]
            
        # Code-phase difference
        elif comb_name == "code_phase":
            try:
                code_phase_1, code_phase_2 = func[comb_name](dset)
            except ValueError:
                log.warn(f"Code-phase difference is not added to dataset. Dual-frequency code and phase observations "
                         f"are needed.")
                continue 
            
            dset.add_float(f"lin.{comb_name}_f1", val=code_phase_1["val"], unit="meter")
            dset.add_float(f"lin.{comb_name}_f2", val=code_phase_2["val"], unit="meter")
            dset.meta.setdefault("linear_combination", dict()).update({f"{comb_name}_f1": code_phase_1["sys_obs"]})
            dset.meta["linear_combination"][f"{comb_name}_f2"] = code_phase_2["sys_obs"]

        # Melbourne-Wuebbena linear combination
        elif comb_name == "melbourne_wuebbena":
            try:
                linear_comb = func[comb_name](dset)
            except ValueError:
                log.warn(f"Melbourne-Wübbena linear combination is not added to dataset. Dual-frequency code and "
                         f"phase observations are needed.")
                continue
            
            dset.add_float(
                f"lin.{comb_name}", 
                val=linear_comb["val"], 
                unit="meter",
            )
            dset.meta.setdefault("linear_combination", dict()).update({f"{comb_name}": linear_comb["sys_obs"]})
        
        else:
            linear_comb = func[comb_name](comb_name, dset)
            for obs_code in linear_comb.keys():
                dset.add_float(
                        f"lin.{comb_name}_{obs_code}", 
                        val=linear_comb[obs_code]["val"], 
                        unit="meter",
                )
                dset.meta.setdefault("linear_combination", dict()).update({f"{comb_name}_{obs_code}": linear_comb[obs_code]["sys_obs"]}) 

            