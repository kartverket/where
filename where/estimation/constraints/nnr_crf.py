import numpy as np

# Midgard imports
from midgard.dev import plugins

# where imports
from where import apriori
from where.lib import config
from where.lib import log


@plugins.register
def nnr_crf(dset, param_names):
    # NNR to CRF

    celestial_reference_frame = config.tech.celestial_reference_frames.list[0]
    crf = apriori.get("crf", time=dset.time, celestial_reference_frame=celestial_reference_frame)
    n = len(param_names)
    H2 = np.zeros((3, n))
    for idx, column in enumerate(param_names):
        if "_src_dir-" not in column:
            continue
        source = column.split("-", maxsplit=1)[-1].split("_")[0]
        if source in crf:
            ra = crf[source].pos.right_ascension
            dec = crf[source].pos.declination
            if column.endswith("_ra"):
                H2[0, idx] = -np.cos(ra) * np.sin(dec) * np.cos(dec)
                H2[1, idx] = -np.sin(ra) * np.sin(dec) * np.cos(dec)
                H2[2, idx] = np.cos(dec) ** 2
            if column.endswith("_dec"):
                H2[0, idx] = np.sin(ra)
                H2[1, idx] = -np.cos(ra)

    if H2.any():
        log.info(f"Applying NNR constraint to {celestial_reference_frame.upper()}")
        # add NNR to CRF constraints
        sigma = np.array([1e-6] * 3)
        return H2, sigma
    else:
        return np.zeros((0, n)), np.zeros(0)
