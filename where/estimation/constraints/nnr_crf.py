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


#     n = len(param_names)
#     d = np.zeros((n, 3))
#     stations = set()
#     #todo: config
#     reference_frame = config.tech.reference_frames.list[0]
#     s = 1.5e-11
#
#     trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=reference_frame)
#
#     # thaller2008: eq 2.51 (skipping scale factor)
#     for idx, column in enumerate(param_names):
#         if "_site_pos-" not in column:
#             continue
#         station = column.split("-", maxsplit=1)[-1].split("_")[0]
#         site_id = dset.meta[station]["site_id"]
#         if site_id in trf:
#             x0, y0, z0 = trf[site_id].pos.trs
#             if column.endswith("_x"):
#                 d[idx, :] = np.array([0, z0, -y0])
#             if column.endswith("_y"):
#                 d[idx, :] = np.array([-z0, 0, x0])
#             if column.endswith("_z"):
#                 d[idx, :] = np.array([y0, -x0, 0])
#             stations.add(station)
#
#     constraint = __name__
#     log.info(f"Applying {constraint} with {', '.join(stations)} from {reference_frame.upper()}")
#     # thaller2008: eq 2.57
#     try:
#         h = np.linalg.inv(d.T @ d) @ d.T
#     except np.linalg.LinAlgError:
#         h = np.zeros((3, n))
#         log.warn(f"Applying {constraint} failed")
#     sigma = np.array([s] * 3)
#     return h, sigma
