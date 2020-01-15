"""Least square estimator

Description:
------------

"""
# Standard library imports
from typing import Dict, List, Union

# External library imports
import numpy as np
from numpy.lib.recfunctions import merge_arrays

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.estimation.estimators._lsq import LsqEstimator
from where.lib import config
from where.lib import log


@plugins.register_named("partial_config_keys")
def partial_config_keys():
    """List the types of partials needed by the estimator

    The LSQ estimator uses only constant parameters. -> TODO: Is that the case?

    Returns:
        Tuple: Strings with names of config keys listing which partial models to run.
    """
    # TODO:
    return ("estimate_constant", "estimate_stochastic")


@plugins.register
def estimate_lsq(dset: "Dataset", partial_vectors: Dict[str, List[str]], obs_noise: np.ndarray) -> None:
    """Estimate with least square method

    The dataset will be updated with least square estimation solution. In addition a convergence status is returned
    via dataset meta variable, which says if the convergence limit is fulfilled.

    Args:
        dset (Dataset):          A dataset containing the data.
        partial_vectors (Dict):  Names and values of the partial derivatives for each partial config key.
        TODO: Used for generation of weight matrix in Kalman filtering? obs_noise (Array):       Observation noise, numpy array with one float value for each observation. 
    """
    # Initialize variables
    estimate_convergence = (
        dset.estimate_convergence if "estimate_convergence" in dset.fields else np.zeros(dset.num_obs)
    )
    convergence_limit = config.tech.convergence_limit.float

    # Organize partial derivatives (parameter vectors) into a matrix
    num_unknowns = len(partial_vectors["estimate_constant"])
    H = np.zeros((dset.num_obs, num_unknowns, 1))
    param_names = list()

    # Constant parameters are simply copied from the partial fields
    for idx, name in enumerate(partial_vectors["estimate_constant"]):
        H[:, idx, 0] = dset["partial." + name][:]
        param_names.append(name)

    # Initialize variables
    z = dset.observed - dset.calc
    W = _weight_matrix(dset)
    x0 = _apriori_values_for_unknowns(dset, param_names)

    # Epochwise estimation or over whole time period
    if config.tech.estimate_epochwise.bool:

        weight_idx_col = 0
        estimate = np.array([])
        for epoch in sorted(set(dset.time.gps.mjd)):
            idx = dset.time.gps.mjd == epoch
            num_obs = dset.time.gps.mjd[idx].size

            # +TODO: Does not work so far. Updating of 'estimate' variable has to be handled.
            ## Skip estimate if convergence limit is fulfilled
            # if np.all(estimate_convergence[idx]):
            #    weight_idx_col+=num_obs
            #    continue
            # -TODO

            # Initialize the estimator for each epoch
            lsq = LsqEstimator(
                H=H[idx],
                z=z[idx],
                x0=x0[idx][0],
                W=W[idx, weight_idx_col : weight_idx_col + num_obs],
                param_names=param_names,
            )

            # Run the estimator
            lsq.estimate()

            # Check if estimated correction is less than convergence limit
            if np.all(lsq.dx <= convergence_limit):
                estimate_convergence[idx] = True

            # Set the starting index for weight matrix
            weight_idx_col += num_obs

            # Append current epoch solution to final estimation solution
            estimate = np.concatenate((estimate, _get_estimate(lsq)), axis=0) if estimate.size else _get_estimate(lsq)
    else:

        # Initialize and run the estimator
        lsq = LsqEstimator(H=H, z=z, x0=x0[0], W=W, param_names=param_names)
        lsq.estimate()

        # Check if estimated correction is less than convergence limit
        if np.all(lsq.dx <= convergence_limit):
            estimate_convergence = np.repeat(True, lsq.num_obs)

        estimate = _get_estimate(lsq)

    log.info(
        f"Variance factor = {np.mean(estimate['estimate_variance_factor']):.4f}, "
        f"degrees of freedom = {np.mean(estimate['estimate_degree_of_freedom']):.1f}"
    )

    # Update the dataset with results from the estimator
    _update_dataset(dset, param_names, estimate, estimate_convergence)


def _get_estimate(lsq: "LsqEstimator") -> np.ndarray:
    """Get estimation solution

    Args:
        lsq:        Least square estimator object.

    Returns:
        Estimation solution.
    """
    # Define field and data types of estimation solution
    dtype = list()
    for section in ["estimate_apriori", "estimate", "estimate_sigma"]:
        dtype.extend([tuple((f"{section}_{param_name}", float)) for param_name in lsq.param_names])
    dtype.extend(
        [
            ("estimate_number_of_observations", float),
            ("estimate_number_of_unknowns", float),
            ("estimate_degree_of_freedom", float),
            ("estimate_variance_factor", float),
            ("residual", float),
        ]
    )

    # Save defined estimation fields in numpy array
    values = [
        tuple(row)
        for row in np.vstack(
            (
                np.repeat(lsq.x0[None, :], lsq.num_obs, axis=0).T,
                np.repeat(lsq.x_hat[None, :], lsq.num_obs, axis=0).T,
                np.repeat(lsq.sigmax[None, :], lsq.num_obs, axis=0).T,
                np.repeat(lsq.num_obs, lsq.num_obs).T,
                np.repeat(lsq.num_unknowns, lsq.num_obs).T,
                np.repeat(lsq.degree_of_freedom, lsq.num_obs).T,
                np.repeat(lsq.sigma0, lsq.num_obs).T,
                lsq.v,
            )
        ).T
    ]

    estimate = np.array(values, dtype=dtype)

    # TODO: Find a more general solution for saving covariance matrix
    if "gnss_site_pos-x" in lsq.param_names:
        covariance_site_pos = _get_gnss_site_pos_covariance(lsq)
        estimate = merge_arrays([estimate, covariance_site_pos], flatten=True)

    if "gnss_rcv_clock" in lsq.param_names:
        variance_rcv_clock = _get_gnss_rcv_clock_variance(lsq)
        estimate = merge_arrays([estimate, variance_rcv_clock], flatten=True)

    return estimate


def _get_gnss_site_pos_covariance(lsq: "LsqEstimator") -> np.ndarray:
    """Get GNSS site position covariance matrix

    Args:
        lsq:  Least square estimator object.

    Returns:
        Covariance matrix of site position.
    """
    # Define field and data types of site position covariance solution
    dtype = [
        ("estimate_cov_site_pos_xx", float),
        ("estimate_cov_site_pos_xy", float),
        ("estimate_cov_site_pos_xz", float),
        ("estimate_cov_site_pos_yy", float),
        ("estimate_cov_site_pos_yz", float),
        ("estimate_cov_site_pos_zz", float),
    ]

    # Get starting index of site position parameters in covariance matrix
    idx = lsq.param_names.index("gnss_site_pos-x")

    # Save defined covariance fields in numpy array
    values = [
        tuple(row)
        for row in np.vstack(
            (
                np.repeat(lsq.Cx[idx, idx], lsq.num_obs, axis=0).T,
                np.repeat(lsq.Cx[idx + 1, idx], lsq.num_obs, axis=0).T,
                np.repeat(lsq.Cx[idx + 2, idx], lsq.num_obs, axis=0).T,
                np.repeat(lsq.Cx[idx + 1, idx + 1], lsq.num_obs, axis=0).T,
                np.repeat(lsq.Cx[idx + 2, idx + 1], lsq.num_obs, axis=0).T,
                np.repeat(lsq.Cx[idx + 2, idx + 2], lsq.num_obs, axis=0).T,
            )
        ).T
    ]

    return np.array(values, dtype=dtype)


def _get_gnss_rcv_clock_variance(lsq: "LsqEstimator") -> np.ndarray:
    """Get GNSS receiver clock variance 

    Args:
        lsq:  Least square estimator object.

    Returns:
        Variance of receiver clock.
    """
    # Define field and data types of receiver clock variance solution
    dtype = [("estimate_cov_rcv_clock_tt", float)]

    # Get starting index of receiver clock parameter in covariance matrix
    idx = lsq.param_names.index("gnss_rcv_clock")

    # Save defined covariance fields in numpy array
    values = [tuple(row) for row in np.vstack((np.repeat(lsq.Cx[idx, idx], lsq.num_obs, axis=0).T,)).T]

    return np.array(values, dtype=dtype)


def _apriori_values_for_unknowns(dset: "Dataset", param_names: List[str]) -> np.ndarray:
    """Get apriori values for unknowns

    If estimates of the unknowns exist, then these values are used as apriori information. Otherwise specified
    apriori values are used for each estimated parameter. 

    Args:
        dset:           A dataset containing the data.
        param_names:    Parameter names.

    Returns:
        Apriori values for unknowns
    """

    x0 = np.zeros((dset.num_obs, len(param_names)))
    apriori_val = {
        "gnss_rcv_clock": np.repeat(0.0, dset.num_obs),
        "gnss_site_pos-x": np.repeat(dset.site_pos.trs.x[0], dset.num_obs),
        "gnss_site_pos-y": np.repeat(dset.site_pos.trs.y[0], dset.num_obs),
        "gnss_site_pos-z": np.repeat(dset.site_pos.trs.z[0], dset.num_obs),
    }

    for idx, param_name in enumerate(param_names):
        field = f"estimate_{param_name}"
        x0[:, idx] = dset[field] if field in dset.fields else apriori_val[param_name]

    return x0


def _weight_matrix(dset: "Dataset") -> np.ndarray:
    """Determine observation weight matrix

    TODO: Distinguish between weighting of code and phase observations

    Args:
        dset: A dataset containing the data.
    """
    sigma0 = config.tech.observation_weight.float
    elevation_weighting = config.tech.elevation_weighting.str

    if elevation_weighting == "sin":
        w_diag = np.sin(dset.site_pos.elevation) / sigma0

    elif elevation_weighting == "sqrtsin":
        w_diag = np.sqrt(np.sin(dset.site_pos.elevation)) / sigma0

    else:
        w_diag = np.ones((dset.num_obs)) / sigma0

    return np.diag(w_diag)


def _update_dataset(
    dset: "Dataset", param_names: List[str], estimate: np.ndarray, estimate_convergence: np.ndarray
) -> None:
    """Update the given dataset with results from the estimation

    Args:
        dset:                   A dataset containing the data.
        param_names:            Parameter names.
        estimate:               Estimation solution.
        estimate_convergence:   Information about estimation convergence for each epoch
    """

    # Update dataset meta
    dset.meta["param_names"] = param_names
    dset.meta["estimate_convergence_status"] = True if np.all(estimate_convergence) else False

    # Update dataset with estimation fields and calculate new residuals
    _add_fields(dset, estimate)
    dset.residual[:] = estimate["residual"].T
    if "estimate_convergence" in dset.fields:
        dset.estimate_convergence[:] = estimate_convergence
    else:
        dset.add_float("estimate_convergence", val=estimate_convergence)

    # Update site position in Dataset
    if "estimate_gnss_site_pos-x" in estimate.dtype.names:
        site_pos = np.array(
            [
                estimate["estimate_gnss_site_pos-x"],
                estimate["estimate_gnss_site_pos-y"],
                estimate["estimate_gnss_site_pos-z"],
            ]
        ).T
        dset.site_pos[:] = site_pos
        # MURKS: Is that needed again? -> dset.site_pos.other(dset.sat_posvel)


def _add_fields(dset: "Dataset", estimate: np.ndarray) -> None:
    """Add fields to the given dataset

    Adds fields for each parameter. 

    Args:
        dset:        A dataset containing the data.
        estimate:    Estimation solution.
    """
    for field in estimate.dtype.names:
        if field in dset.fields:
            dset[field][:] = estimate[field]
        else:
            dset.add_float(
                field, val=estimate[field], unit="meter", write_level="operational"  # TODO: Better solution for unit?
            )
