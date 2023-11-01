"""General Kalman filter

Description:
------------

Kalman filter and the modified Bryson-Frazier smoother is covered in the book "Factorization Methods for Discrete
Sequential Estimation" by Bierman :cite:`bierman2006`. However, there are some typos but a corrected version of the
algorithm is listed in :cite:`gibbs2011`.

"""

# External library imports
import numpy as np
import h5py

# Standard library import
import os

# Midgard imports
from midgard.math.unit import Unit

# Where imports
from where.lib import log
from where.lib import config


class KalmanFilter(object):
    """A general Kalman filter

    See for instance https://en.wikipedia.org/wiki/Kalman_filter#Details for information about a general Kalman
    filter. We use the Modified Bryson-Frazier smoother, which is described at
    https://en.wikipedia.org/wiki/Kalman_filter#Modified_Bryson.E2.80.93Frazier_smoother

    Notation:

    h:                 Partial derivatives                                          # num_obs x n x 1
    x:                 Predicted state estimate (x-tilde)                           # num_obs x n x 1
    x_hat:             Updated state estimate (x-hat)                               # num_obs x n x 1
    sigma:             Residual covariance                                          # num_obs
    innovation:        Measurement residual                                         # num_obs
    z:                 Observed residual                                            # num_obs
    r:                 Observation noise covariance                                 # num_obs
    Q:                 Process noise covariance                                     # dict()
    k:                 Kalman gain                                                  # num_obs x n x 1
    p:                 Predicted estimate covariance (p-tilde)                      # n x n  (not stored)
    p_hat:             Updated estimate covariance (p-hat)                          # num_obs x n x n
    phi:               State transition                                             # n x n  (not stored)

    x_smooth:          Smoothed state estimates                                     # num_obs x n x 1
    lam:               (lambda)                                                     # num_obs x n x 1
    """

    def __init__(self, h, z=None, apriori_stdev=None, phi=None, r=None, Q=None, param_names=None):
        """Initialize the Kalman filter

        Args:
            h (Numpy array):               Partial derivatives          (num_obs x n x 1)
            z (Numpy array):               Observations                 (num_obs)
            apriori_stdev (Numpy array):   Apriori standard deviation   (n)
            phi (Numpy array):             State transition             (num_obs x n x n)
            r (Numpy array):               Observation noise covariance (num_obs)
            Q (Numpy array):               Process noise covariance     (num_obs x n x n)
        """
        self.h = h
        self.num_obs, self.n, _ = self.h.shape
        self.apriori_stdev = np.ones(self.n) if apriori_stdev is None else apriori_stdev

        self.z = np.zeros((self.num_obs)) if z is None else z
        self.phi = np.eye(self.n).repeat(self.num_obs).reshape(self.n, self.n, -1).T if phi is None else phi
        self.r = np.ones((self.num_obs)) if r is None else r
        self.Q = dict() if Q is None else Q
        self.x_hat = np.zeros((self.num_obs, self.n, 1))
        self.x_hat_ferr = np.zeros((self.num_obs, self.n))
        self.x_smooth = np.zeros((self.num_obs, self.n, 1))
        self.param_names = param_names if param_names else []
        self.innovation = np.zeros(self.num_obs)
        self.sigma = np.zeros(self.num_obs)

        self.p_hat_file_path = config.files.path("output_covariance_matrix")
        self.p_hat_file = h5py.File(self.p_hat_file_path, "w")
        self.p_hat_file.attrs["labels"] = ", ".join(self.param_names)
        self.p_hat_file.close()

    def filter(self):
        """Run the Kalman filter forward and backward
        """
        # Initialize
        x_tilde = np.zeros((self.n, 1))
        p_tilde = np.diag(self.apriori_stdev ** 2)
        k = np.zeros((self.num_obs, self.n, 1))
        lam = np.zeros((self.n, 1))

        # Makes calculations easier to read (and gives a slight speed-up)
        h = self.h
        z = self.z
        phi = self.phi
        r = self.r
        Q = self.Q
        x_hat = self.x_hat
        x_smooth = self.x_smooth
        I = np.eye(self.n)
        innovation = self.innovation
        sigma = self.sigma

        # Run filter forward over all observations
        for epoch in range(self.num_obs):
            innovation[epoch] = z[epoch] - h[epoch].T @ x_tilde
            sigma[epoch] = (h[epoch].T @ p_tilde @ h[epoch]) + r[epoch]
            k[epoch] = p_tilde @ h[epoch] / sigma[epoch]
            x_hat[epoch] = x_tilde + k[epoch] * innovation[epoch]
            p_hat = (I - k[epoch] @ h[epoch].T) @ p_tilde

            if isinstance(phi[epoch], int):
                # phi is identity matrix. Save computation time by skipping multiplication with identity matrix
                x_tilde = x_hat[epoch]
                p_tilde = p_hat
            else:
                x_tilde = phi[epoch] @ x_hat[epoch]
                p_tilde = phi[epoch] @ p_hat @ phi[epoch].T

            for (idx1, idx2), noise in Q.get(epoch, {}).items():
                p_tilde[idx1, idx2] += noise

            self._set_p_hat(epoch, p_hat)
            self.x_hat_ferr[epoch, :] = np.sqrt(np.diagonal(p_hat))

        # Run smoother backwards over all observations
        for epoch in range(self.num_obs - 1, -1, -1):
            # TODO smooth covariance matrix
            p_hat = self._get_p_hat(epoch)
            x_smooth[epoch] = x_hat[epoch] + p_hat.T @ lam
            if isinstance(phi[epoch - 1], int):
                # phi is identity matrix. Remove it from equation to speed up computation
                lam = (
                    h[epoch] * innovation[epoch] / sigma[epoch]
                    + (I - k[epoch] @ h[epoch].T).T @ lam
                )
            else:
                lam = (
                    phi[epoch - 1].T @ h[epoch] * innovation[epoch] / sigma[epoch]
                    + phi[epoch - 1].T @ (I - k[epoch] @ h[epoch].T).T @ lam
                )

    def update_dataset(self, dset, param_names, normal_idx, num_unknowns):
        """Update the given dataset with results from the filtering

        Args:
            dset (Dataset):       The dataset.
            param_names (List):   Strings with names of parameters. Used to form field names.
            normal_idx (Slice):   Slice denoting which parameters should be used for the normal equations.
            num_unknowns (Int):   Number of unknowns.
        """
        # Update dataset with state and estimation fields and calculate new residuals
        self._add_fields(dset, param_names)
        dset.residual[:] = dset.est - (dset.obs - dset.calc)
        num_unknowns += dset.meta.get("num_clock_coeff", 0)

        # Calculate normal equations, and add statistics about estimation to dataset
        N, b = self._normal_equations(normal_idx, dset.num_obs - 1)
        g = self.x_hat[dset.num_obs - 1, normal_idx, :]
        deg_freedom = dset.num_obs - num_unknowns
        v = dset.residual[:, None]
        P = np.diag(1 / self.r[: dset.num_obs])
        sq_sum_residuals = (v.T @ P @ v).item()
        sq_sum_omc_terms = (2 * b.T @ g - g.T @ N @ g).item()
        variance_factor = sq_sum_residuals / deg_freedom if deg_freedom != 0 else np.inf
        log.info(f"Variance factor = {variance_factor:.4f}, degrees of freedom = {deg_freedom:d}")

        # Report and set analysis status if there are too few degrees of freedom
        if deg_freedom < 1:
            log.error(f"Degrees of freedom is {deg_freedom} < 1. Estimate fewer parameters")
            if dset.meta.get("analysis_status") == "unchecked":
                dset.meta["analysis_status"] = "too few degrees of freedom"

        else:
            if dset.meta.get("analysis_status") == "too few degrees of freedom":
                dset.meta["analysis_status"] = "unchecked"

        # Report and set analysis status if there are too few stations
        # TODO: if vlbi_site_pos in state_vector and num_stations < 3
        estimate_site_pos = np.char.startswith(np.array(param_names, dtype=str), "vlbi_site_pos").any()
        if len(dset.unique("station")) < 3 and estimate_site_pos:
            log.warn(f"Too few stations {len(dset.unique('station'))} < 3. Do not estimate station positions.")
            # if dset.meta.get("analysis_status") == "unchecked":
            # dset.meta["analysis_status"] = "needs custom state vector"
        elif len(dset.unique("station")) < 3 and estimate_site_pos:
            if dset.meta.get("analysis_status") == "needs custom state vector":
                dset.meta["analysis_status"] = "unchecked"
        # Update config
        cfg_vars = dset.vars.copy()
        cfg_vars.pop("rundate")
        with config.update_tech_config(dset.analysis["rundate"], cfg_vars.pop("pipeline"), **cfg_vars) as cfg:
            cfg.update("analysis_status", "status", dset.meta.get("analysis_status", ""), source=__file__)

        # Add information to dset.meta
        dset.meta.add("number of observations", dset.num_obs, section="statistics")
        dset.meta.add("number of unknowns", num_unknowns, section="statistics")
        dset.meta.add("square sum of residuals", sq_sum_residuals, section="statistics")
        dset.meta.add("degrees of freedom", deg_freedom, section="statistics")
        dset.meta.add("variance factor", variance_factor, section="statistics")
        dset.meta.add("weighted square sum of o-c", sq_sum_residuals + sq_sum_omc_terms, section="statistics")
        dset.meta.add("matrix", N.tolist(), section="normal equation")
        dset.meta.add("vector", b[:, 0].tolist(), section="normal equation")
        dset.meta.add("names", param_names[normal_idx], section="normal equation")
        dset.meta.add(
            "unit", [config.tech[f.split("-")[0]].unit.str for f in param_names[normal_idx]], section="normal equation"
        )

    def cleanup(self):
        if not config.tech.keep_covariance_file.bool:
            os.remove(self.p_hat_file_path)

    def _add_fields(self, dset, param_names):
        """Add fields to the given dataset

        Adds fields for state vectors and estimate vectors for each parameter. Parameters with names ending with an
        underscore, `_`, are not added to the dataset.

        Args:
            dset (Dataset):       The dataset.
            param_names (List):   Strings with names of parameters. Used to form field names.

        """
        # Delete values from previous iterations
        if "state" in dset.fields:
            del dset.state

        if "estimate" in dset.fields:
            del dset.estimate

        for idx, param_name in enumerate(param_names):
            if param_name.endswith("_"):
                continue

            # State vectors
            fieldname = f"state.{param_name}"
            fieldname_sigma = fieldname + "_sigma"
            value = self.x_smooth[: dset.num_obs, idx, 0]
            value_sigma = np.sqrt(self.x_hat_ferr[: dset.num_obs, idx])

            # Convert values to the display unit. It corresponds to "meter per <unit of partial>"
            partial_unit = dset.unit("partial.{}".format(param_name))
            to_unit = dset.meta["display_units"][param_name]
            from_unit = f"meter/({partial_unit[0]})"
            factor = Unit(from_unit, to_unit)
            dset.meta.add(param_name, factor, section="display_factors")
            dset.add_float(fieldname, val=value * factor, unit=to_unit, write_level="operational")

            # Convert values to the display unit. It corresponds to "meter per <unit of partial>"
            partial_unit = dset.unit("partial.{}".format(param_name))
            to_unit = dset.meta["display_units"][param_name]
            from_unit = f"meter/({partial_unit[0]})"
            factor = Unit(from_unit, to_unit)
            dset.meta.add(param_name, factor, section="display_factors")
            dset.add_float(fieldname_sigma, val=value_sigma * factor, unit=to_unit, write_level="operational")

            # Estimate vectors
            fieldname = f"estimate.{param_name}"
            value = self.h[: dset.num_obs, idx, 0] * self.x_smooth[: dset.num_obs, idx, 0]
            dset.add_float(fieldname, val=value, unit="meter", write_level="analysis")

        value = (self.x_smooth.transpose(0, 2, 1) @ self.h)[: dset.num_obs, 0, 0]
        fieldname = "est"
        if fieldname in dset.fields:
            dset[fieldname][:] = value
        else:
            dset.add_float(fieldname, val=value, unit="meter", write_level="operational")

    def _normal_equations(self, normal_idx, last_obs):
        """Calculate normal equations corresponding to the filter results

        Args:
            normal_idx (Slice):  A slice denoting which columns should be used for the normal equations.

        Returns:
            Tuple of Numpy arrays: Normal matrix (n x n) and Normal vector (n x 1).
        """
        p_hat_last = self._get_p_hat(last_obs)
        p_tilde_0 = np.diag(self.apriori_stdev ** 2)

        pg = p_hat_last[normal_idx, normal_idx]
        cg = p_tilde_0[normal_idx, normal_idx]
        g = self.x_hat[last_obs, normal_idx, :]

        pg_inv = np.linalg.inv(pg)
        cg_inv = np.linalg.inv(cg)

        N = pg_inv - cg_inv
        b = pg_inv @ g

        # test
        if False:
            stat_idx = slice(normal_idx.stop, self.n, None)
            R = np.diag(self.r[: last_obs + 1])
            H_L = self.h[: last_obs + 1, stat_idx, 0]
            c_L = p_tilde_0[stat_idx, stat_idx]
            R_tilde = H_L @ c_L @ H_L.T + R
            R_tilde_inv = np.linalg.inv(R_tilde)

            H_g = self.h[: last_obs + 1, normal_idx, 0]

            NN = H_g.T @ R_tilde_inv @ H_g
            bb = H_g.T @ R_tilde_inv @ self.z[: last_obs + 1]

            print("Normal matrix allclose:", np.allclose(NN, N), "Max diff:", np.max(np.abs(NN - N)))
            print("Normal vector allclose:", np.allclose(bb, b[:, 0]), "Max diff:", np.max(np.abs(bb - b[:, 0])))
            # import IPython; IPython.embed()
        # end test
        return N, b

    def _set_p_hat(self, epoch, data):
        with h5py.File(self.p_hat_file_path, "a") as fid:
            fid.create_dataset(str(epoch), data=data)

    def _get_p_hat(self, epoch):
        with h5py.File(self.p_hat_file_path, "r") as fid:
            return fid[str(epoch)][...]
