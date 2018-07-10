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

# Where imports
from where.lib.unit import unit
from where.lib import log
from where.lib import config
from where.lib import files


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
        self.p_hat_file_path = files.path("output_covariance_matrix")

        with h5py.File(self.p_hat_file_path, "w") as p_hat_file:
            p_hat_file.attrs["labels"] = ", ".join(self.param_names)

    def filter(self):
        """Run the Kalman filter forward and backward
        """
        # Initialize
        x_tilde = np.zeros((self.n, 1))
        p_tilde = np.diag(self.apriori_stdev ** 2)
        sigma = np.zeros(self.num_obs)
        innovation = np.zeros(self.num_obs)
        k = np.zeros((self.num_obs, self.n, 1))
        lam = np.zeros((self.n, 1))

        # Makes calculations easier to read (and gives a slight speed-up)
        h = self.h
        z = self.z
        phi = self.phi
        r = self.r
        Q = self.Q
        x_hat = self.x_hat
        x_hat_ferr = self.x_hat_ferr
        x_smooth = self.x_smooth
        I = np.eye(self.n)

        # Run filter forward over all observations
        for epoch in range(self.num_obs):
            innovation[epoch] = z[epoch] - h[epoch].T @ x_tilde
            sigma[epoch] = (h[epoch].T @ p_tilde @ h[epoch]) + r[epoch]
            k[epoch] = p_tilde @ h[epoch] / sigma[epoch]
            x_hat[epoch] = x_tilde + k[epoch] * innovation[epoch]
            p_hat = (I - k[epoch] @ h[epoch].T) @ p_tilde

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
            lam = (
                phi[epoch - 1].T @ h[epoch]
                * innovation[epoch]
                / sigma[epoch]
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
        dset.residual[:] = dset.obs + dset.estimate - dset.calc
        num_unknowns += dset.meta["num_clock_coeff"] if "num_clock_coeff" in dset.meta else 0

        # Calculate normal equations, and add statistics about estimation to dataset
        N, b = self._normal_equations(normal_idx, dset.num_obs - 1)
        g = self.x_hat[dset.num_obs - 1, normal_idx, :]
        deg_freedom = dset.num_obs - num_unknowns
        v = dset.residual[:, None]
        P = np.diag(1 / self.r[:dset.num_obs])
        sq_sum_residuals = np.asscalar(v.T @ P @ v)
        sq_sum_omc_terms = np.asscalar(2 * b.T @ g - g.T @ N @ g)
        variance_factor = sq_sum_residuals / deg_freedom
        log.info("Variance factor = {:.4f}", variance_factor)

        # Report and set analysis status if there are too few degrees of freedom
        if deg_freedom < 1:
            log.error(f"Degrees of freedom is {deg_freedom} < 1. Estimate fewer parameters")
            if dset.meta["analysis_status"] == "unchecked":
                dset.meta["analysis_status"] = "too few degrees of freedom"

                # Update config
                # with config.update_tech_config(dset.rundate, dset.vars["tech"], dset.vars["session"]) as cfg:
                # cfg.update("analysis_status", "status", dset.meta["analysis_status"], source=__file__)
        else:
            if dset.meta["analysis_status"] == "too few degrees of freedom":
                dset.meta["analysis_status"] = "unchecked"

                # Update config
                # with config.update_tech_config(dset.rundate, dset.vars["tech"], dset.vars["session"]) as cfg:
                # cfg.update("analysis_status", "status", dset.meta["analysis_status"], source=__file__)

        # Report and set analysis status if there are too few stations
        # TODO: if vlbi_site_pos in state_vector and num_stations < 3
        if len(dset.unique("station")) < 3:
            log.error(f"Two few stations {len(dset.unique('station'))} < 3. Do not estimate station positions.")
            if dset.meta["analysis_status"] == "unchecked":
                dset.meta["analysis_status"] = "needs custom state vector"

        print(dset.meta["analysis_status"])
        # Update config
        with config.update_tech_config(dset.rundate, dset.vars["tech"], dset.vars["session"]) as cfg:
            cfg.update("analysis_status", "status", dset.meta["analysis_status"], source=__file__)

        # Add information to dset.meta
        dset.add_to_meta("statistics", "number of observations", dset.num_obs)
        dset.add_to_meta("statistics", "number of unknowns", num_unknowns)
        dset.add_to_meta("statistics", "square sum of residuals", sq_sum_residuals)
        dset.add_to_meta("statistics", "degrees of freedom", deg_freedom)
        dset.add_to_meta("statistics", "variance factor", variance_factor)
        dset.add_to_meta("statistics", "weighted square sum of o-c", sq_sum_residuals + sq_sum_omc_terms)
        dset.add_to_meta("normal equation", "matrix", N.tolist())
        dset.add_to_meta("normal equation", "vector", b[:, 0].tolist())
        dset.add_to_meta("normal equation", "names", param_names[normal_idx])
        dset.add_to_meta(
            "normal equation", "unit", [config.tech[f.split("-")[0]].unit.str for f in param_names[normal_idx]]
        )

        # TODO should this be here?
        log.info("Solving normal equations")
        names = dset.meta["normal equation"]["names"]
        n = len(names)
        d = np.zeros((n, 6))
        stations = set()
        reference_frame = config.tech.reference_frames.list[0]

        from where import apriori

        trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=reference_frame)
        celestial_reference_frame = config.tech.celestial_reference_frames.list[0]
        crf = apriori.get("crf", celestial_reference_frames=celestial_reference_frame, session=dset.dataset_name)

        # thaller2008: eq 2.51 (skipping scale factor)
        for idx, column in enumerate(names):
            if "_site_pos-" not in column:
                continue
            station = column.split("-", maxsplit=1)[-1].split("_")[0]
            site_id = dset.meta[station]["site_id"]
            if site_id in trf:
                x0, y0, z0 = trf[site_id].pos.itrs  # TODO: Take units into account
                if column.endswith("_x"):
                    d[idx, :] = np.array([1, 0, 0, 0, z0, -y0])
                if column.endswith("_y"):
                    d[idx, :] = np.array([0, 1, 0, -z0, 0, x0])
                if column.endswith("_z"):
                    d[idx, :] = np.array([0, 0, 1, y0, -x0, 0])
                stations.add(station)

        log.info("Applying NNT/NNR with {} from {}", ", ".join(stations), reference_frame.upper())
        # thaller2008: eq 2.57
        try:
            H = np.linalg.inv(d.T @ d) @ d.T
        except np.linalg.LinAlgError:
            H = np.zeros((6, n))

        sigmas = [0.0001] * 3 + [1.5e-11] * 3

        # NNR to CRF
        H2 = np.zeros((3, n))
        for idx, column in enumerate(names):
            if "_src_dir-" not in column:
                continue
            source = column.split("-", maxsplit=1)[-1].split("_")[0]
            if source in crf:
                ra = crf[source].pos.crs[0]
                dec = crf[source].pos.crs[1]
                if column.endswith("_ra"):
                    H2[0, idx] = -np.cos(ra) * np.sin(dec) * np.cos(dec)
                    H2[1, idx] = -np.sin(ra) * np.sin(dec) * np.cos(dec)
                    H2[2, idx] = np.cos(dec) ** 2
                if column.endswith("_dec"):
                    H2[0, idx] = np.sin(ra)
                    H2[1, idx] = -np.cos(ra)

        if H2.any():
            log.info("Applying NNR constraint to {}", celestial_reference_frame.upper())
            # add NNR to CRF constraints
            H = np.concatenate((H, H2))
            sigmas = sigmas + [1e-6] * 3
        # thaller2008: eq 2.45
        P_h = np.diag(1 / np.array(sigmas) ** 2)

        # thaller2008: eq 2.58
        N_h = N + H.T @ P_h @ H

        # solve neq
        N_h_inv = np.linalg.inv(N_h)
        x = N_h_inv @ b

        # Covariance: thaller2008: eq 2.16
        Q_xx = variance_factor ** 2 * N_h_inv

        dset.add_to_meta("normal equation", "solution", x[:, 0].tolist())
        dset.add_to_meta("normal equation", "covariance", Q_xx.tolist())

    def cleanup(self):
        if not config.tech.keep_covariance_file.bool:
            files.delete_file(self.p_hat_file_path)

    def _add_fields(self, dset, param_names):
        """Add fields to the given dataset

        Adds fields for state vectors and estimate vectors for each parameter. Parameters with names ending with an
        underscore, `_`, are not added to the dataset.

        Args:
            dset (Dataset):       The dataset.
            param_names (List):   Strings with names of parameters. Used to form field names.

        """
        for idx, param_name in enumerate(param_names):
            if param_name.endswith("_"):
                continue

            # State vectors
            fieldname = "{}_{}".format("state", param_name)
            fieldname_sigma = fieldname + "_sigma"
            value = self.x_smooth[:dset.num_obs, idx, 0]
            value_sigma = np.sqrt(self.x_hat_ferr[:dset.num_obs, idx])

            if fieldname in dset.fields:
                dset[fieldname][:] = value * dset.meta["display_factors"][param_name]
            else:
                # Convert values to the display unit. It corresponds to "meter per <unit of partial>"
                partial_unit = dset.unit("partial_{}".format(param_name))
                display_unit = dset.meta["display_units"][param_name]
                factor = unit("meter / ({})".format(partial_unit), display_unit)
                dset.add_to_meta("display_factors", param_name, factor)
                dset.add_float(
                    fieldname, table="state", val=value * factor, unit=display_unit, write_level="operational"
                )

            if fieldname_sigma in dset.fields:
                dset[fieldname_sigma][:] = value * dset.meta["display_factors"][param_name]
            else:
                # Convert values to the display unit. It corresponds to "meter per <unit of partial>"
                partial_unit = dset.unit("partial_{}".format(param_name))
                display_unit = dset.meta["display_units"][param_name]
                factor = unit("meter / ({})".format(partial_unit), display_unit)
                dset.add_to_meta("display_factors", param_name, factor)
                dset.add_float(
                    fieldname_sigma,
                    table="state",
                    val=value_sigma * factor,
                    unit=display_unit,
                    write_level="operational",
                )

            # Estimate vectors
            fieldname = "{}_{}".format("estimate", param_name)
            value = self.h[:dset.num_obs, idx, 0] * self.x_smooth[:dset.num_obs, idx, 0]
            if fieldname in dset.fields:
                dset[fieldname][:] = value
            else:
                dset.add_float(fieldname, table="estimate", val=value, unit="meter", write_level="analysis")

        value = (self.x_smooth.transpose(0, 2, 1) @ self.h)[:dset.num_obs, 0, 0]
        if "estimate" in dset.fields:
            dset.estimate[:] = value
        else:
            dset.add_float("estimate", val=value, unit="meter", write_level="operational")

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
            R = np.diag(self.r[:last_obs + 1])
            H_L = self.h[:last_obs + 1, stat_idx, 0]
            c_L = p_tilde_0[stat_idx, stat_idx]
            R_tilde = H_L @ c_L @ H_L.T + R
            R_tilde_inv = np.linalg.inv(R_tilde)

            H_g = self.h[:last_obs + 1, normal_idx, 0]

            NN = H_g.T @ R_tilde_inv @ H_g
            bb = H_g.T @ R_tilde_inv @ self.z[:last_obs + 1]
            print("Normal matrix all close:", np.allclose(NN, N, atol=0.01), "Max diff:", np.max(np.abs(NN - N)))
            print(
                "Normal vector all close:",
                np.allclose(bb, b[:, 0], atol=0.01),
                "Max diff:",
                np.max(np.abs(bb - b[:, 0])),
            )
            # import IPython; IPython.embed()
        # end test
        return N, b

    def _set_p_hat(self, epoch, data):
        with h5py.File(self.p_hat_file_path) as p_hat_file:
            p_hat_file.create_dataset(str(epoch), data=data)

    def _get_p_hat(self, epoch):
        with h5py.File(self.p_hat_file_path) as p_hat_file:
            return p_hat_file[str(epoch)][...]
