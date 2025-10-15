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
from midgard.math.unit import Unit

# Where imports
from where import apriori
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
    """Solve normal equation matrix stored in dset.meta
    
    The solution is stored in dset.meta 
    
    Reference:
    Daniela Thaller Phd Thesis 2008: "Inter-technique combination based on homogeneous normal equation 
    systems including station coordinates, Earth orientation and troposphere parameters"
    """
    log.info("Solving normal equations")
    names = dset.meta["normal equation"]["names"]
    n = len(names)
    abs_param_weight = np.zeros(n)
    rel_param_weight = np.zeros(n)
    H = np.zeros((6, n))
    H_sigma = np.ones(6)

    reference_frame = config.tech.minimum_trf.reference_frame.str or config.tech.reference_frames.list[0]    

    B, stations = _compute_helmert_matrix(dset, reference_frame)
    
    # NNT/NNR to TRF
    if "minimum_trf" in config.tech.neq_constraints.list:
        if len(stations) >= 3:
            H, H_sigma, constraints = _apply_minimum_trf(H, H_sigma, B, stations)
            log.info(f"Applying {'/'.join(constraints)} with {', '.join(stations)} from {reference_frame.upper()}")
        else:
            log.info(
                f"Too few stations to use minimum contraints from {reference_frame.upper()}. Using absolute constraints for station positions."
                )
            # Automatic absolute constraints: Too few stations to use NNT/NNR
            for idx, column in enumerate(names):
                if "_site_pos-" not in column:
                    continue
                parameter, _,subparameter = column.partition("-")
                station = subparameter.rsplit("_", maxsplit=1)[0]
                weight = config.tech[parameter].neq_abs_constraint_weight.float
                abs_param_weight[idx] += 1 / weight ** 2  # 1/meters**2
        

    # NNR to CRF
    if "minimum_crf" in config.tech.neq_constraints.list:
        H, H_sigma = _apply_minimum_crf(dset, H, H_sigma)

    # thaller2008: eq 2.45
    P_h = np.diag(1 / np.array(H_sigma) ** 2)

    # Free network constraints: thaller2008: eq 2.58
    N = np.array(dset.meta["normal equation"]["matrix"])
    N_h = N + H.T @ P_h @ H

    # Automatic absolute constraints; Baselines with too few obs?
    for idx, column in enumerate(names):
        if "_baseline-" not in column:
            continue
        parameter, _,subparameter = column.partition("-")
        baseline = subparameter.rsplit("_", maxsplit=1)[0]
        if dset.num(baseline=baseline) < 5:
            weight = config.tech[parameter].neq_abs_constraint_weight.float
            abs_param_weight[idx] = 1 / weight ** 2  # 1/meters**2
            log.info(f"Too few observations for baseline {baseline}. Constrained to a priori value")
            continue

    # Check for configured absolute constraints
    for idx, column in enumerate(names):
        parameter = column.split("-", maxsplit=1)[0]
        if not "neq_abs_constraint" in config.tech[parameter]:
            continue
        
        if config.tech[parameter].neq_abs_constraint.bool:
            log.info(f"Applying absolute constraints for {parameter}")
            weight = config.tech[parameter].neq_abs_constraint_weight.float
            abs_param_weight[idx] = 1 / weight**2

    # Apply absolute constraints. thaller2008: eq.2.49
    N_h += np.diag(abs_param_weight)

    # Check for configured relative constraints
    part = ""
    num_params = 0
    previous_parameter = ""
    H_rel_const = np.zeros((n,n))
    for idx, column in enumerate(names):
        parameter, subparameter = column.split("-", maxsplit=1)
        
        if config.tech[parameter].neq_rate_constraint.bool and "neq_rate_constraint" in config.tech[parameter]:
            
            weight = config.tech[parameter].neq_rate_constraint_weight.float
            unit = config.tech[parameter].unit.str
            knot_interval = config.tech[parameter].knot_interval.int * Unit.seconds2hours
            weight_per_interval = weight * knot_interval
            rel_param_weight[idx] = 1 /weight_per_interval**2
            new_part, epoch = subparameter.split(":", maxsplit=1)
            num_params += 1
            if new_part != part and part != "":
                log.info(f"Applying relative constraints of {weight:6.4f}{unit}/h for {previous_parameter} {part}")
                H_part = np.diag(np.ones(num_params)) - np.diag(np.ones(num_params-1), 1)
                # Skip the last row in H_part
                H_rel_const[idx-(num_params-1):idx, (idx-(num_params-1)):idx+1] = H_part[:-1, :]
                rel_param_weight[idx] = 0
                num_params = 0    
            part = new_part
        else: 
            if num_params != 0 and part != "":
                # Apply the final relative constraint when there is a new parameter in the list
                log.info(f"Applying relative constraints of {weight:6.4f}{unit}/h for {previous_parameter} {part}")
                H_part = np.diag(np.ones(num_params)) - np.diag(np.ones(num_params-1), 1)
                # Skip the last row in H_part
                H_rel_const[idx-(num_params-1):idx, (idx-(num_params-1)):idx+1] = H_part[:-1, :]
                rel_param_weight[idx] = 0
                num_params = 0
                part = ""
        previous_parameter = parameter

    if num_params != 0 and part != "":
        # Apply the final relative constraint when this was the last parameter in the list
        log.info(f"Applying relative constraints of {weight:6.4f}{unit}/h for {previous_parameter} {part}")
        H_part = np.diag(np.ones(num_params)) - np.diag(np.ones(num_params-1), 1)
        # Skip the last row in H_part
        H_rel_const[idx-(num_params-1):idx, (idx-(num_params-1)):idx+1] = H_part[:-1, :]
        rel_param_weight[idx] = 0
        num_params = 0
        part = ""
                
    # Apply relative constraints. thaller2008: eq.2.60
    N_h += H_rel_const.T @ np.diag(rel_param_weight) @ H_rel_const
        
    # solve neq
    N_h_inv = np.linalg.inv(N_h)
    b = np.array(dset.meta["normal equation"]["vector"])
    x = N_h_inv @ b[:, None]

    num_abs_constraints = np.sum(abs_param_weight != 0)
    num_rel_constraints = np.sum(rel_param_weight != 0)
    
    # Update statistics for solution after constraints are added
    v_c = H @ x  # Sinex Format description appendix equation 10
    dset.meta["statistics"]["square sum of residuals"] += (v_c.T @ P_h @ v_c).item()
    dset.meta["statistics"]["degrees of freedom"] += len(H) + int(num_abs_constraints) + int(num_rel_constraints)
    dset.meta["statistics"]["variance factor"] = dset.meta["statistics"]["square sum of residuals"] / dset.meta["statistics"]["degrees of freedom"]

    variance_factor = dset.meta["statistics"]["variance factor"]
    deg_freedom = dset.meta["statistics"]["degrees of freedom"]
    log.info(f"Updated variance factor = {variance_factor:.4f}, degrees of freedom = {deg_freedom:d}")

    # Covariance: thaller2008: eq 2.16
    Q_xx = variance_factor ** 2 * N_h_inv

    dset.meta.add("solution", x[:, 0].tolist(), section="normal equation")
    dset.meta.add("covariance", Q_xx.tolist(), section="normal equation")

    _compute_helmert_parameters(dset, B, stations, "neq")


def _compute_helmert_matrix(dset, reference_frame):
    """Compute the Jacobian matrix for the 7-helmert parameters

    The matrix is computed with the following order of the parameters:
        T_x, T_y, T_z, D, R_x, R_y, R_z
    
    Args:
        dset (Dataset):               Dataset with model data
        reference_frame (String):     Reference frame for the transformation parameters

    Returns:
        B (array):    Jacobian matrix for the Helmert parameters
        stations:     List of stations used to form the matrix
    """
    names = dset.meta["normal equation"]["names"]
    n = len(names)
    B = np.zeros((n, 7))
    stations = list()
    
    trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=reference_frame)
    skip_stations = config.tech.minimum_trf.skip_stations.list
    # thaller2008: eq 2.51
    for idx, column in enumerate(names):
        if "_site_pos-" not in column:
            continue
        station = column.split("-", maxsplit=1)[-1].rsplit("_", maxsplit=1)[0]
        site_id = dset.meta["station"][station]["site_id"]
        if station not in skip_stations:
            try:
                x0, y0, z0 = trf[site_id].pos.trs
            except KeyError:
                # Station is not defined in the given reference frame for the given time
                continue
            #IERS 2010 Conventions eq. 4.8
            if column.endswith("_x"):
                B[idx, :] = np.array([1, 0, 0, x0, 0, z0, -y0])
                stations.append(station)
            if column.endswith("_y"):
                B[idx, :] = np.array([0, 1, 0, y0, -z0, 0, x0])
            if column.endswith("_z"):
                B[idx, :] = np.array([0, 0, 1, z0, y0, -x0, 0])
    return B, stations 
           
def _apply_minimum_trf(H, H_sigma, B, stations):
    """Compute minimum constraints to the TRF
    
    Args:
        H (array):       Initial Jacobian matrix for the constraint
        H_sigma (array): Initial weights for the constraint
        B (array):       Jacobian matrix for the 7 Helmert parameters
        stations (Set):  List of station names used for the Helmert parameters
        
    Returns:
        H (array):       Updated Jacobian matrix for the constraint
        H_sigma (array): Updated weights for the constraint
    """
    constraints = []
    nns = config.tech.minimum_trf["nns"].bool
    nnt = config.tech.minimum_trf["nnt"].bool
    nnr = config.tech.minimum_trf["nnr"].bool
    
    trf_nns_unit = config.tech.minimum_trf.nns_unit.str
    trf_nns_sigma = config.tech.minimum_trf.nns_sigma.float * Unit(trf_nns_unit, "unit") # Convert to unit

    trf_nnt_unit = config.tech.minimum_trf.nnt_unit.str 
    trf_nnt_sigma = config.tech.minimum_trf.nnt_sigma.float * Unit(trf_nnt_unit, "meter") # Convert to meter
    
    trf_nnr_unit = config.tech.minimum_trf.nnr_unit.str 
    trf_nnr_sigma = config.tech.minimum_trf.nnr_sigma.float * Unit(trf_nnr_unit, "rad") # Convert to radians

    if not nns:
        # Remove scale factor
        d = np.delete(B, 3, axis=1)
    else:
        d = B

    try:
        # thaller2008: eq 2.57
        H = np.linalg.inv(d.T @ d) @ d.T
    except np.linalg.LinAlgError:
        log.warn(f"Unable to invert matrix for NNR/NNT constraints")

    if nnt and nnr and not nns:
        H_sigma[0:3] = trf_nnt_sigma
        H_sigma[3:6] = trf_nnr_sigma
        constraints.append("NNT")
        constraints.append("NNR")
    elif nnt and nnr and nns:
        H_sigma[0:3] = trf_nnt_sigma
        H_sigma[3:6] = trf_nnr_sigma
        H_sigma = np.insert(H_sigma, 3, trf_nns_sigma)
        constraints.append("NNT")
        constraints.append("NNR")
        constraints.append("NNS")
    elif nnt and not nnr and not nns:
        # This scenario is not tested
        H = d[:, 0:3]
        H_sigma[0:3] = trf_nnt_sigma
        H_sigma = H_sigma[0:3]
        constraints.append("NNT")
    elif nnt and not nnr and nns:
        # This scenario is not tested
        H = d[:, 0:4]
        H_sigma[0:3] = trf_nnt_sigma
        H_sigma = H_sigma[0:3]
        H_sigma = np.insert(H_sigma, 3, trf_nns_sigma)
        constraints.append("NNT")
        constraints.append("NNS")
    elif not nnt and nnr and not nns:
        # This scenario is not tested
        H = d[:, 3:6]
        H_sigma[3:6] = trf_nnt_sigma
        H_sigma = H_sigma[3:6]
        constraints.append("NNR")
    elif not nnt and nnr and nns:
        # This scenario is not tested
        H = d[:, 3:7]
        H_sigma[3:6] = trf_nnt_sigma
        H_sigma = H_sigma[3:6]
        H_sigma = np.insert(H_sigma, 0, trf_nns_sigma)
        constraints.append("NNR")
        constraints.append("NNS")
    elif not nnt and not nnr and not nns:
        # This scenario is not tested
        #H = np.zeros((n, 0))
        pass
    elif not nnt and not nnr and nns:
        # This scenario is not tested
        H = d[3, 3]
        H_sigma = np.array([trf_nns_sigma])
        constraints.append("NNS")
    else:
        #H = np.zeros((n, 0))
        log.warn(f"Unknown constraints. Not applying.")
    
    return H, H_sigma, constraints
        
def _apply_minimum_crf(dset, H, H_sigma):
    """Compute minimum constraints to the CRF
    
    Args:
        dset (Dataset):  Dataset with model data
        H (array):       Initial Jacobian matrix for the constraint
        H_sigma (array): Initial weights for the constraint
        
    Returns:
        H (array):       Updated Jacobian matrix for the constraint
        H_sigma (array): Updated weights for the constraint
    """
    names = dset.meta["normal equation"]["names"]
    n = len(names)
    
    frame = config.tech.minimum_crf.reference_frame.str or config.tech.celestial_reference_frames.list[0]
    crf = apriori.get("crf", time=dset.time, celestial_reference_frames=frame)
    skip_sources = config.tech.minimum_crf.skip_sources.list # TODO: only defining sources
    H2 = np.zeros((3, len(names)))
    for idx, column in enumerate(names):
        if "_src_dir-" not in column:
            continue
        parameter, _,subparameter = column.partition("-")
        source = subparameter.split("_")[0]
        # Source names are saved internally in Where with the letters "dot" instead of the character "." 
        # which is originally in the source name for some sources
        source = source.replace("dot", ".")
        if source in crf and source not in skip_sources:
            ra = crf[source].pos.right_ascension
            dec = crf[source].pos.declination
            if dset.num(source=source) < 5:
                weight = config.tech[parameter].neq_abs_constraint_weight.float
                abs_param_weight[idx] = 1 / weight ** 2  # 1/radians**2
                if column.endswith("_ra"):
                    # Avoid logging the same message twice for both ra and dec
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
        log.info(f"Applying NNR constraint to {frame.upper()}")
        # add NNR to CRF constraints
        H = np.concatenate((H, H2))
        crf_nnr_unit = config.tech.minimum_crf.unit.str 
        crf_nnr_sigma = config.tech.minimum_crf.sigma.float * Unit(crf_nnr_unit, "rad") # Convert to radians
        H_sigma2 = np.full(3, fill_value=crf_nnr_sigma)
        H_sigma = np.concatenate((H_sigma, H_sigma2))

    return H, H_sigma

def _compute_helmert_parameters(dset, B, stations, source=None):
    """ Estimate Helmert parameters for solution with regards to apriori reference frame
    
    The estimated parameters are added to dset.meta
    """
    def _get_estimated_pos(station, param, source=None):
        """Get estimated position either from solved NEQ or Kalman filter solution"""
        name = f"{config.tech.master_section.name}_site_pos-{station}_{param}"
        if source == "neq":
            idx = dset.meta["normal equation"]["names"].index(name)
            return dset.meta["normal equation"]["solution"][idx]
        else:
            return np.mean(dset.state[name])

    reference_frame = config.tech.reference_frames.list[0]
    trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=reference_frame)
    pos_ref = []
    pos_est = []
    for station in stations:
        site_id = dset.meta["station"][station]["site_id"]
        x0, y0, z0 = trf[site_id].pos.trs
        pos_ref.append(x0)
        pos_ref.append(y0)
        pos_ref.append(z0)
        x1 = x0 + _get_estimated_pos(station, "x", source)
        y1 = y0 + _get_estimated_pos(station, "y", source)
        z1 = z0 + _get_estimated_pos(station, "z", source)
        pos_est.append(x1)
        pos_est.append(y1)
        pos_est.append(z1)

    keep_idx = np.sum(B, axis=1) != 0
    B = B[keep_idx]
    try:
        hp = np.linalg.inv(B.T @ B) @ B.T @ (np.array(pos_est) - np.array(pos_ref))
    except np.linalg.LinAlgError:
        hp = np.full(len(B.T), fill_value=np.nan)
        
    section = f"{source}_helmert" if source else "helmert"
    fields = ["T_X", "T_Y", "T_Z", "scale", "alpha", "beta", "gamma"]
    width = max([len(f) for f in fields])
    factors = [Unit.m2mm, Unit.m2mm, Unit.m2mm, Unit.unit2ppb, Unit.rad2mas, Unit.rad2mas, Unit.rad2mas]
    units = ["mm", "mm", "mm", "ppb", "mas", "mas", "mas"]
    if source:
        log.info(f"Session {source.upper()} Helmert parameters with regards to {reference_frame.upper()}")
    else:
        log.info(f"Session Kalman Filter Helmert parameters with regards to {reference_frame.upper()}")
    for i, (field, factor, unit) in enumerate(zip(fields, factors, units)):
        dset.meta.add(field, (float(hp[i] * factor), unit), section=section)
        log.info(
            f"{field:{width}} = {dset.meta[section][field][0]: 6.4f} [{dset.meta[section][field][1]}]"
        )
