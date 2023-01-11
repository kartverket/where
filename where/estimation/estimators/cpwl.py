"""Continuous PieceWise Linear estimator

Description:
------------

"""

# External library imports
import numpy as np
import scipy.sparse

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.estimation.estimators._kalman import KalmanFilter
from where.lib import config
from where.lib import log


@plugins.register_named("partial_config_keys")
def partial_config_keys():
    """List the types of partials needed by the estimator

    The CPWL estimator uses both constant and stochastic parameters.

    Returns:
        Tuple: Strings with names of config keys listing which partial models to run.
    """
    return ("estimate_constant", "estimate_stochastic")


@plugins.register
def estimate_cpwl(dset, partial_vectors, obs_noise):
    """Estimate with continuous piecewise linear functions

    TODO: Describe phi and Q

    Args:
        dset (Dataset):          Model run data.
        partial_vectors (Dict):  Names and values of the partial derivatives for each partial config key.
        obs_noise (Array):       Observation noise, numpy array with one float value for each observation.
    """
    # Organize partial derivatives (state vectors) into a matrix
    n_constant = len(partial_vectors["estimate_constant"])
    n_stochastic = len(partial_vectors["estimate_stochastic"])
    n = n_constant + 2 * n_stochastic
    num_unknowns = n
    num_obs = dset.num_obs
    h = np.zeros((num_obs, n, 1))
    param_names = list()

    # Constant parameters are simply copied from the partial fields
    for idx, name in enumerate(partial_vectors["estimate_constant"]):
        h[:, idx, 0] = dset["partial." + name][:]
        param_names.append(name)

    # Stochastic parameters are estimated as CPWL functions by adding a rate parameter
    for idx, name in enumerate(partial_vectors["estimate_stochastic"]):
        h[:, n_constant + idx * 2, 0] = dset["partial." + name][:]
        param_names.extend([name, name + "_rate_"])  # Trailing underscore in rate_ means field is not added to dset

    # Read information about parameters from config files
    ref_time = np.ones(n) * dset.time.utc[0].mjd
    knot_interval = np.ones(n) * np.inf
    process_noise = np.zeros(n)
    apriori_stdev = np.empty(n)

    constant_params = {c.split("-")[0] for c in partial_vectors["estimate_constant"]}
    for param in constant_params:
        idx = np.array([c.startswith(param + "-") for c in param_names])
        apriori_stdev[idx] = config.tech[param].apriori_stdev.float

    stochastic_params = {c.split("-")[0] for c in partial_vectors["estimate_stochastic"]}
    for param in stochastic_params:
        # Set default knot_interval
        intervals = config.tech[param].knot_interval.list
        const_idx = np.array([c.startswith(param + "-") for c in param_names])
        rate_idx = np.array([c.startswith(param + "-") and c.endswith("rate_") for c in param_names])

        knot_interval[rate_idx] = float(intervals.pop(0)) * Unit.seconds2day
        for interval in intervals:
            # (Potentially) overwrite with station specific knot_interval
            sta, _, seconds = interval.partition(":")
            rate_idx_sta = np.array(
                [c.startswith(param + "-") and c.endswith("rate_") and sta in c for c in param_names]
            )
            knot_interval[rate_idx_sta] = float(seconds) * Unit.seconds2day
        process_noise[rate_idx] = config.tech[param].process_noise.float

        apriori_stdev[const_idx] = config.tech[param].apriori_stdev.float
        apriori_stdev[rate_idx] = config.tech[param].apriori_rate_stdev.float  # Rate parameters

    # Initialize variables
    z = dset.obs - dset.calc
    # phi = np.repeat(np.eye(n)[None, :, :], num_obs, axis=0)
    phi = list()

    delta_phi = np.eye(n, k=1)
    delta_phi[:, :n_constant] = 0
    delta_phi[:, n_constant::2] = 0

    Q = dict()

    for epoch in range(num_obs - 1):
        # TODO: Check that 24 is correct here (and use unit instead)
        delta_t = (dset.time.utc[epoch + 1].mjd - dset.time.utc[epoch].mjd) * 24
        # phi[epoch] += delta_phi * delta_t
        phi.append(scipy.sparse.csr_matrix(np.eye(n) + delta_phi * delta_t))

        idx = np.logical_and(process_noise, dset.time.utc[epoch + 1].mjd > ref_time + knot_interval)
        indicies = np.where(idx)[0]
        Q[epoch] = {(i, i): process_noise[i] ** 2 for i in indicies}

        ref_time[idx] += knot_interval[idx] * ((dset.time.utc[epoch + 1].mjd - ref_time[idx]) // knot_interval[idx])
        num_unknowns += int(sum(idx))

    phi.append(scipy.sparse.csr_matrix(np.eye(n)))

    # Add pseudo-observations
    constraints = config.tech.get(key="estimate_constraint", default="").as_list(split_re=", *")
    if "minimum_trf" in constraints:
        frame = config.tech.minimum_trf.reference_frame.str or config.tech.reference_frames.list[0]
        trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=frame)
        d = np.zeros((n, 6))
        stations = set()
        skip_stations = config.tech.minimum_trf.skip_stations.list

        for idx, column in enumerate(param_names):
            if "_site_pos-" not in column:
                continue
            station = column.split("-", maxsplit=1)[-1].rsplit("_", maxsplit=1)[0]
            key = dset.meta["station"][station]["site_id"]
            if key in trf and station not in skip_stations:
                x0, y0, z0 = trf[key].pos.trs  # TODO: Take units into account
                if column.endswith("_x"):
                    d[idx, :] = np.array([1, 0, 0, 0, z0, -y0])
                if column.endswith("_y"):
                    d[idx, :] = np.array([0, 1, 0, -z0, 0, x0])
                if column.endswith("_z"):
                    d[idx, :] = np.array([0, 0, 1, y0, -x0, 0])
                stations.add(station)

        # TODO deal with slr_site_pos etc
        if "vlbi_site_pos" in constant_params:
            trf_constraints = [c.upper() for c in ["nnt", "nnr"] if config.tech.minimum_trf[c].bool]
            log.info(
                f"Applying {'/'.join(trf_constraints)} with {', '.join(stations)} from {frame.upper()}"
            )
            nnt = config.tech.minimum_trf["nnt"].bool
            nnr = config.tech.minimum_trf["nnr"].bool
            
            nnt_unit = config.tech.minimum_trf.nnt_unit.str 
            nnt_sigma = config.tech.minimum_trf.nnt_sigma.float * Unit(nnt_unit, "meter") # Convert to meter
            
            nnr_unit = config.tech.minimum_trf.nnr_unit.str 
            nnr_sigma = config.tech.minimum_trf.nnr_sigma.float * Unit(nnr_unit, "rad") # Convert to radians
            
            if nnt and nnr:
                obs_noise = np.hstack((obs_noise, np.array([nnt_sigma ** 2] * 3 + [nnr_sigma ** 2] * 3))).T
            elif nnt and not nnr:
                d = d[:, 0:3]
                obs_noise = np.hstack((obs_noise, np.array([nnt_sigma ** 2] * 3))).T
            elif not nnt and nnr:
                d = d[:, 3:6]
                obs_noise = np.hstack((obs_noise, np.array([nnr_sigma ** 2] * 3))).T
            elif not nnt and not nnr:
                d = np.zeros((n, 0))
            else:
                d = np.zeros((n, 0))
                log.warn(f"Unknown constraints: {'/'.join(constraints).upper()}. Not applying.")

        num_constraints = d.shape[1]
        try:
            h = np.vstack((h, (np.linalg.inv(d.T @ d) @ d.T)[:, :, None]))
        except np.linalg.linalg.LinAlgError:
            pass

        if "vlbi_src_dir" in constant_params:
            if "minimum_crf" in constraints:
                frame = config.tech.minimum_crf.reference_frame.str or config.tech.celestial_reference_frames.list[0]
                crf = apriori.get("crf", time=dset.time, celestial_reference_frames=frame)
                skip_sources = config.tech.minimum_crf.skip_sources.list # TODO: only defining sources
                # NNR to CRF
                log.info(f"Applying NNR constraint to {frame.upper()}")
                H2 = np.zeros((3, n))
                for idx, column in enumerate(param_names):
                    if "_src_dir-" not in column:
                        continue
                    source = column.split("-", maxsplit=1)[-1].split("_")[0]
                    if source in crf and source not in skip_sources:
                        ra = crf[source].pos.right_ascension
                        dec = crf[source].pos.declination
                        if column.endswith("_ra"):
                            H2[0, idx] = -np.cos(ra) * np.sin(dec) * np.cos(dec)
                            H2[1, idx] = -np.sin(ra) * np.sin(dec) * np.cos(dec)
                            H2[2, idx] = np.cos(dec) ** 2
                        if column.endswith("_dec"):
                            H2[0, idx] = np.sin(ra)
                            H2[1, idx] = -np.cos(ra)
    
                nnr_unit = config.tech.minimum_crf.unit.str 
                nnr_sigma = config.tech.minimum_crf.sigma.float * Unit(nnr_unit, "rad") # Convert to radians
                obs_noise = np.hstack((obs_noise, np.array([nnr_sigma ** 2] * 3)))
                num_constraints += 3
                h = np.vstack((h, H2[:, :, None]))

        z = np.hstack((z, np.zeros(num_constraints))).T
        # phi = np.vstack((phi, np.repeat(np.eye(n)[None, :, :], num_constraints, axis=0)))
        phi = phi + [scipy.sparse.csr_matrix(np.eye(n))] * num_constraints

    # Initialize and run the Kalman filter
    kalman = KalmanFilter(h, z=z, apriori_stdev=apriori_stdev, phi=phi, r=obs_noise, Q=Q, param_names=param_names)
    kalman.filter()

    # Update the dataset with results from the filter
    kalman.update_dataset(dset, param_names=param_names, normal_idx=slice(0, n_constant), num_unknowns=num_unknowns)
    kalman.cleanup()
