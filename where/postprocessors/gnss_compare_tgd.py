"""Compare total/broadcast group delay given from navigation messages against other differential code bias (DCB) sources

Description:
------------
Information to calculate the difference between total/broadcast group delay (TGD/BGD) and DCBs taken from: 
Wang, N. et al (2019): "Quality assessment of GPS, Galileo and BeiDou-2/3 satellite broadcast group delays"

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
from where.lib import exceptions, log

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
    """Compare total/broadcast group delay (TGD/BGD) given from navigation messages against other differential code 
    bias (DCB) sources

    Depending on the given GNSS navigation data following fields can be added to Dataset:

      | Field                | Type          | Unit | Description                                                     |
      | :------------------- | :------------ | :--- | :-------------------------------------------------------------- |
      | {solution}_dcb       | numpy.ndarray | s    | BGD/TGD calculated based on postprocessed DCBs                  |
      | {solution}_dcb_mean  | numpy.ndarray | s    | BGD/TGD calculated based on postprocessed DCBs. In addition the |
      |                      |               |      | mean over the data period is subtracted from the solution (zero |
      |                      |               |      | mean difference) for each GNSS seperately.                      |
      | {solution}_mean      | numpy.ndarray | s    | The mean over the data period is subtracted from the BGD/TDG    |
      |                      |               |      | solution (zero mean difference) for each GNSS seperately.       |
      | {solution}_diff      | numpy.ndarray | s    | BGD/TGD compared to given DCBs                                  |
      | {solution}_diff_mean | numpy.ndarray | s    | BGD/TGD compared to given DCBs. In addition the mean over the   |
      |                      |               |      | data period is subtracted from the solution (zero mean          |
      |                      |               |      | difference) for each GNSS seperately.                           |

    whereby following {solution} are possible:

      | Solution   | Description                         |
      | :----------| :---------------------------------- |
      | bgd_e1_e5a | Galileo E1-E5a BGD                  | 
      | bgd_e1_e5b | Galileo E1-E5a BGD                  | 
      | tgd        | Total group delay (e.g. from GPS)   | 
      | tgd_b1_b3  | BeiDou B1-B3 TGD                    |
      | tgd_b2_b3  | BeiDou B2-B3 TGD                    |
        
    Args:
        dset:     A Dataset containing model data.
    """
    
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
            fields_to_add = [f"{field}_mean", f"{field}_dcb", f"{field}_dcb_mean", f"{field}_diff", f"{field}_diff_mean"]
            if f"{field}_dcb" not in dset.fields:
                tmp_nan = np.empty(dset.num_obs)
                tmp_nan[:] = np.nan
                log.info(f"Add TGD/BGD comparison and DCB fields to dataset: {field}_dcb, {field}_diff, {field}_diff_mean")    
                for field_to_add in fields_to_add:
                    dset.add_float(field_to_add, val=tmp_nan.copy(), unit="second")

            if tgd_def[sys][field].field  not in dset.fields:
                dset.add_float(tgd_def[sys][field].field, val=tmp_nan.copy(), unit="second")

            # Add values to added dataset fields
            for sat in dset.unique("satellite", idx=idx_sys):
                idx = dset.filter(satellite=sat)
       
                try:
                    dcb_val = dcb.get_dsb(sat, tgd_def[sys][field].dcb, dset.analysis["rundate"])["estimate"]
                except exceptions.MissingDataError as error:
                    # Note: DCBs are not availbale in bias file.
                    log.debug(error)
                    continue
                    
                dset[tgd_def[sys][field].field][idx] = dcb_val
                
                if sys in ["E", "G", "I", "J"]:
                    f1 = getattr(enums, "gnss_freq_" + sys)[tgd_def[sys][field].f1].value
                    f2 = getattr(enums, "gnss_freq_" + sys)[tgd_def[sys][field].f2].value
                    dset[f"{field}_dcb"][idx] = dcb_val/(1-(f1**2/f2**2))
                    dset[f"{field}_diff"][idx] = dset[field][idx] - dset[f"{field}_dcb"][idx]
                    
                elif sys == "C":
                    dset[f"{field}_diff"][idx] = dset[field][idx] - dcb_val

            # Check if DCBs for given field has been available
            if np.all(np.isnan(dset[f"{field}_dcb"])):
                log.warn(f"No DCBs available for {enums.gnss_id_to_name[sys].value}. Related dataset '{field}' "
                         f"fields will be deleted.")
                for field_to_add in fields_to_add:
                    del dset[field_to_add]
                continue

            # Determine zero mean difference between TGD and post-processed DCBs
            dcb_sys_mean = np.mean(dset[f"{field}_dcb"][idx_sys])
            tgd_sys_mean = np.mean(dset[field][idx_sys])
            dset[f"{field}_dcb_mean"][idx_sys] = dset[f"{field}_dcb"][idx_sys] - dcb_sys_mean
            dset[f"{field}_mean"][idx_sys] = dset[field][idx_sys] - tgd_sys_mean
            dset[f"{field}_diff_mean"][idx_sys] = dset[f"{field}_diff"][idx_sys] - tgd_sys_mean + dcb_sys_mean
    
