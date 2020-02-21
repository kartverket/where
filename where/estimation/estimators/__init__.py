"""Framework for running estimators

Description:
------------

Each estimator should be defined in a separate .py-file. The function inside the .py-file that
should be called need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def estimate_least_square(dset):
        ...

The decorated function will be called with a single parameter, ``dset`` which contains a
:class:`~where.data.dataset.Dataset` with data that can be used in the estimation. If the estimator needs the partial
derivatives it should obtain them by itself calling the :func:`where.estimation.partial_vectors` function::

    from where import estimation
    partial_vectors = estimation.partial_vectors(dset, 'estimate', 'estimate_stochastic')

"""
# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log


def call(config_key, dset, partial_vectors, obs_noise):
    """Call an estimator

    Args:
        config_key (String):     Config key specifying the name of the estimator.
        dset (Dataset):          Model run data.
        partial_vectors (Dict):  Names and values of the partial derivatives for each partial config key.
        obs_noise (Array):       Observation noise, numpy array with one float value for each observation.
    """
    estimator_name = config.tech[config_key].str
    if estimator_name:
        plugins.call(
            package_name=__name__,
            plugin_name=estimator_name,
            dset=dset,
            partial_vectors=partial_vectors,
            obs_noise=obs_noise,
        )


def partial_config_keys(estimator_config_key):
    """Find which partials that a given estimator requires.

    Finds the partial config keys by calling the function registered as 'partial_config_keys' on the given estimator.

    Args:
        estimator_config_key (String):   Config key specifying the name of the estimator.

    Returns:
        Tuple: Strings with names of config keys listing which partial models to run.
    """
    estimator_name = config.tech[estimator_config_key].str
    return plugins.call(package_name=__name__, plugin_name=estimator_name, part="partial_config_keys")


def solve_neq(dset):
    log.info("Solving normal equations")
    names = dset.meta["normal equation"]["names"]
    n = len(names)
    d = np.zeros((n, 6))
    fix_param_weight = np.zeros(n)
    H = np.zeros((6, n))
    stations = set()
    from where import apriori

    reference_frame = config.tech.reference_frames.list[0]
    trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=reference_frame)
    # thaller2008: eq 2.51 (skipping scale factor)
    for idx, column in enumerate(names):
        if "_site_pos-" not in column:
            continue
        station = column.split("-", maxsplit=1)[-1].rsplit("_", maxsplit=1)[0]
        site_id = dset.meta[station]["site_id"]
        if site_id in trf:
            x0, y0, z0 = trf[site_id].pos.trs
            if column.endswith("_x"):
                d[idx, :] = np.array([1, 0, 0, 0, z0, -y0])
            if column.endswith("_y"):
                d[idx, :] = np.array([0, 1, 0, -z0, 0, x0])
            if column.endswith("_z"):
                d[idx, :] = np.array([0, 0, 1, y0, -x0, 0])
            stations.add(station)

    if len(stations) >= 3:
        try:
            # thaller2008: eq 2.57
            H = np.linalg.inv(d.T @ d) @ d.T
            log.info(f"Applying NNT/NNR with {', '.join(stations)} from {reference_frame.upper()}")
        except np.linalg.LinAlgError:
            log.warn(f"Unable to invert matrix for NNR/NNT constraints")
    else:
        log.info(
            f"Too few stations to use NNR/NNT contraints from {reference_frame.upper()}. Using absolute constraints for station positions."
        )
        # Too few stations to use NNT/NNR?
        for idx, column in enumerate(names):
            if "_site_pos-" not in column:
                continue
            station = column.split("-", maxsplit=1)[-1].rsplit("_", maxsplit=1)[0]
            fix_param_weight[idx] = 1 / (1e-6) ** 2  # 1/meters**2

    sigmas = [0.0001] * 3 + [1.5e-11] * 3

    # NNR to CRF
    if "celestial_reference_frames" in config.tech.master_section:
        celestial_reference_frame = config.tech.celestial_reference_frames.list[0]
        crf = apriori.get("crf", time=dset.time, celestial_reference_frames=celestial_reference_frame)
        H2 = np.zeros((3, n))
        for idx, column in enumerate(names):
            if "_src_dir-" not in column:
                continue
            source = column.split("-", maxsplit=1)[-1].split("_")[0]
            if source in crf:
                ra = crf[source].pos.right_ascension
                dec = crf[source].pos.declination
                if dset.num(source=source) < 5:
                    fix_param_weight[idx] = 1 / (1e-12) ** 2  # 1/radians**2
                    if column.endswith("_ra"):
                        log.info(
                            f"Too few observations for source {source}. Using absolute constraints for source positions."
                        )
                    continue

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
            H = np.concatenate((H, H2))
            sigmas = sigmas + [1e-6] * 3

    # thaller2008: eq 2.45
    P_h = np.diag(1 / np.array(sigmas) ** 2)

    # Free network constraints: thaller2008: eq 2.58
    N = np.array(dset.meta["normal equation"]["matrix"])
    N_h = N + H.T @ P_h @ H

    # Baselines with too few obs?
    for idx, column in enumerate(names):
        if "_baseline-" not in column:
            continue
        baseline = column.split("-", maxsplit=1)[-1].rsplit("_", maxsplit=1)[0]
        if dset.num(baseline=baseline) < 5:
            fix_param_weight[idx] = 1 / (1e-6) ** 2  # 1/meters**2
            log.info(f"Too few observations for baseline {baseline}. Constrained to a priori value")
            continue

    # Absolute constraints (on sources with too few observations): thaller2008: eq.2.49
    N_h += np.diag(fix_param_weight)

    # solve neq
    N_h_inv = np.linalg.inv(N_h)
    b = np.array(dset.meta["normal equation"]["vector"])
    x = N_h_inv @ b[:, None]

    # Covariance: thaller2008: eq 2.16
    variance_factor = dset.meta["statistics"]["variance factor"]
    Q_xx = variance_factor ** 2 * N_h_inv

    dset.meta.add("solution", x[:, 0].tolist(), section="normal equation")
    dset.meta.add("covariance", Q_xx.tolist(), section="normal equation")
