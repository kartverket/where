"""Compare total/broadcast group delay given from navigation messages against other differential code bias (DCB) sources

Description:
------------
Information to calculate the difference between total/broadcast group delay (TGD/BGD) and DCBs taken from: 
Wang, N. et al (2019): "Quality assessment of GPS, Galileo and BeiDou-2/3 satellite broadcast group delays"

Depending on the given GNSS navigation data following fields can be added to Dataset:

| Field          | Type          | Unit | Description                                              |
|----------------|---------------|------|----------------------------------------------------------|
| bgd_e1_e5a_diff | numpy.ndarray | s    | Galileo E1-E5a BGD compared to given DCBs                | 
| bgd_e1_e5b_diff | numpy.ndarray | s    | Galileo E1-E5a BGD compared to given DCBs                | 
| tgd_diff        | numpy.ndarray | s    | Total group delay (e.g. from GPS) compared to given DCBs | 
| tgd_b1_b3_diff  | numpy.ndarray | s    | BeiDou B1-B3 TGD compared to given DCBs                  |
| tgd_b2_b3_diff  | numpy.ndarray | s    | BeiDou B2-B3 TGD compared to given DCBs                  |

"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.collections import enums
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


BiasField = namedtuple(
    "BiasField", ["f1", "f2", "dcb", "field"]
)
BiasField.__new__.__defaults__ = (None,) * len(BiasField._fields)
BiasField.__doc__ = """A convenience class for defining a necessary parameters to determine bias based on DCBs

    Args:
        f1 (str):         1st frequency name
        f2 (str):         2nd frequency name
        dcb (str):        Differential code bias (DCB)
        field (str):      DCB field name
    """

@plugins.register
def gnss_compare_tgd(dset: "Dataset") -> None:
    
    tgd_def = {
        "C": {
            "tgd_b1_b3": BiasField("B1", "B3", "C2I-C6I", "dcb_c2i_c6i"),
            "tgd_b2_b3": BiasField("B2", "B3", "C7I-C6I", "dcb_c7i_c6i"),
        },
        "E": {
            "bgd_e1_e5a":  BiasField("E1", "E5a", "C1C-C5Q", "dcb_c1c_c5q"),
            "bgd_e1_e5b":  BiasField("E1", "E5b", "C1C-C7Q", "dcb_c1c_c7q"),
        },
        "G": {
            "tgd": BiasField("L1", "L2", "C1W-C2W", "dcb_c1w_c2w"),
        },
        "J": {
            "tgd": BiasField("L1", "L2", "C1X-C2X", "dcb_c1x_c2x"),  #TODO: Is that correct?
        },
        #"I": {
        #    "tgd": BiasField("L1", "L2", "C1W-C2W", "dcb_c1w_c2w"),  
        #},
        #"R": {
        #    "tgd": BiasField("L1", "L2", "C1P-C2P", "dcb_c1p_c2p"),
        #},
        
    }

    dcb = apriori.get("gnss_bias", rundate=dset.analysis["rundate"])
    

    for sys in dset.unique("system"):
        
        # Skip if GNSS is not defined
        if not sys in tgd_def.keys():
            continue
        
        idx_sys = dset.filter(system=sys)
            
        for field in tgd_def[sys].keys():
            
            # Skip if TGD is not defined in Dataset
            if field not in dset.fields:
                continue
                       
            # Initialize dataset fields
            if f"{field}_dcb" not in dset.fields:
                tmp_nan = np.empty(dset.num_obs)
                tmp_nan[:] = np.nan
                log.info(f"Add TGD/BGD comparison and DCB fields to dataset: {field}_dcb, {field}_diff, {field}_diff_mean")
                dset.add_float(f"{field}_mean", val=tmp_nan.copy(), unit="second")
                dset.add_float(f"{field}_dcb", val=tmp_nan.copy(), unit="second")
                dset.add_float(f"{field}_dcb_mean", val=tmp_nan.copy(), unit="second")
                dset.add_float(f"{field}_diff", val=tmp_nan.copy(), unit="second")
                dset.add_float(f"{field}_diff_mean", val=tmp_nan.copy(), unit="second")

            if tgd_def[sys][field].field  not in dset.fields:
                dset.add_float(tgd_def[sys][field].field, val=tmp_nan.copy(), unit="second")

            # Add values to added dataset fields
            for sat in dset.unique("satellite", idx=idx_sys):
                idx = dset.filter(satellite=sat)
       
                dcb_val = dcb.get_dsb(sat, tgd_def[sys][field].dcb, dset.analysis["rundate"])["estimate"]
                dset[tgd_def[sys][field].field][idx] = dcb_val
                
                if sys in ["E", "G", "I", "J"]:
                    f1 = getattr(enums, "gnss_freq_" + sys)[tgd_def[sys][field].f1].value
                    f2 = getattr(enums, "gnss_freq_" + sys)[tgd_def[sys][field].f2].value
                    dset[f"{field}_dcb"][idx] = dcb_val/(1-(f1**2/f2**2))
                    dset[f"{field}_diff"][idx] = dset[field][idx] - dset[f"{field}_dcb"][idx]
                    
                elif sys == "C":
                    dset[f"{field}_diff"][idx] = dset[field][idx] - dcb_val

            # Determine zero mean difference between TGD and post-processed DCBs
            dcb_sys_mean = np.mean(dset[f"{field}_dcb"][idx_sys])
            tgd_sys_mean = np.mean(dset[field][idx_sys])
            dset[f"{field}_dcb_mean"][idx_sys] = dset[f"{field}_dcb"][idx_sys] - dcb_sys_mean
            dset[f"{field}_mean"][idx_sys] = dset[field][idx_sys] - tgd_sys_mean
            dset[f"{field}_diff_mean"][idx_sys] = dset[f"{field}_diff"][idx_sys] - tgd_sys_mean + dcb_sys_mean
    

