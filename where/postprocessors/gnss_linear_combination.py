"""Add GNSS linar observation combinations to dataset

Description:
------------

Depending on the configuration, following linear combination can be added:
    geometry_free
    ionosphere_free
    wide_lane
    narrow_lane
    melbourne_wuebbena
    code_multipath

"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.lib.gnss import linear_combination, linear_combination_cmc, linear_combination_melbourne

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
            
            dset.add_float(f"lin.{comb_name}", val=cmc1["val"], unit="meter")
            dset.add_float(f"lin.{comb_name}", val=cmc2["val"], unit="meter")
        
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
        
        else:
            linear_comb = func[comb_name](comb_name, dset)
            for obs_code in linear_comb.keys():
                dset.add_float(
                        f"lin.{comb_name}", 
                        val=linear_comb[obs_code]["val"], 
                        unit="meter",
                )
 

            