"""General least square estimator

Description:
------------

"""
# Standard library imports
from typing import List, Union

# External library imports
import numpy as np

# Where imports
from where.lib import log


class LsqEstimator(object):
    """A general least square estimator

    Notation:

    H:                 Design matrix with partial derivatives           # num_obs x num_unknowns
    z:                 Observed residual (obs - calc)                   # num_obs x 1
    W:                 Observation weight matrix                        # num_obs x num_obs
    x0:                Apriori values of estimated parameters           # num_unknowns x 1
    x_hat:             Estimated solution (x-hat = x0 + dx)             # num_unknowns x 1
    dx:                Estimated corrections                            # num_unknowns x 1
    v:                 Estimated residuals                              # num_obs x 1
    sigma0:            Estimated standard deviation of unit weight      # 1
    sigmax:            Standard deviation of the unknowns               # num_unknowns x 1
    N:                 Normal equation                                  # num_unknowns x num_unknowns
    Qx:                Cofactor matrix of the unknowns                  # num_unknowns x num_unknowns
    Cx:                Covariance matrix of the unknowns                # num_unknowns x num_unknowns
    Ql:                Cofactor matrix of the estimated observations    # num_obs x num_obs
    Cl:                Covariance matrix of the estimated observations  # num_obs x num_obs
    Qv:                Cofactor matrix of the residuals                 # num_obs x num_obs
    r:                 Redundancy                                       # num_obs x 1
    """

    def __init__(
        self,
        H: np.ndarray,
        z: Union[None, np.ndarray] = None,
        x0: Union[None, np.ndarray] = None,
        W: Union[None, np.ndarray] = None,
        param_names: Union[None, List[str]] = None,
    ) -> None:
        """Initialize the Kalman filter

        Args:
            H:            Design matrix with partial derivatives  (num_obs x num_unknowns)
            z:            Observed residual                       (num_obs)
            x0:           Apriori values of estimated parameters  (num_unknowns x 1)
            W:            Observation weight matrix               (num_obs x num_obs)
            param_names:  Parameter names                         (num_unknowns)
            
        """
        self.H = np.squeeze(H.T, axis=0).T
        self.num_obs, self.num_unknowns = self.H.shape
        self.degree_of_freedom = self.num_obs - self.num_unknowns

        if self.degree_of_freedom < 0:
            log.error(f"Degree of freedom is {self.degree_of_freedom} < 0. Estimate fewer parameters.")

        self.z = np.zeros((self.num_obs)) if z is None else z
        self.x0 = np.zeros((self.num_unknowns)) if x0 is None else x0
        self.W = np.eye(self.num_obs) if W is None else W  # Initialze as identity matrix, if not given
        self.param_names = param_names if param_names else []

        self.dx = np.zeros((self.num_unknowns))
        self.x_hat = np.zeros((self.num_unknowns))
        self.N = np.zeros((self.num_unknowns, self.num_unknowns))
        self.Cx = np.zeros((self.num_unknowns, self.num_unknowns))
        self.Qx = np.zeros((self.num_unknowns, self.num_unknowns))
        self.v = np.zeros((self.num_obs))
        self.sigma0 = None
        self.sigmax = np.zeros((self.num_unknowns))

    def estimate(self) -> None:
        """Run the least square estimator
        """
        # Makes calculations easier to read (and gives a slight speed-up)
        H = self.H
        z = self.z
        W = self.W
        x0 = self.x0

        # Solution of normal equations
        self.N = H.T @ W @ H
        if not np.isfinite(np.linalg.cond(self.N)):
            log.warn("Error by computing the inverse of normal equation matrix N.")
        self.dx = np.linalg.inv(self.N) @ H.T @ W @ z
        self.x_hat = x0 - self.dx

        # Estimated residuals
        self.v = H @ self.dx - z

        # Estimated standard deviation of unit weight
        self.sigma0 = np.sqrt(self.v.T @ W @ self.v / self.degree_of_freedom)

        # Cofactor matrix of the unknowns
        self.Qx = np.linalg.inv(self.N)

        # Covariance matrix of the unknowns
        self.Cx = self.sigma0 ** 2 * self.Qx

        # Standard deviation of the unknowns
        self.sigmax = np.sqrt(np.diag(self.Cx))

        # Estimated observations
        l0 = H @ self.dx

        # Cofactor matrix of the estimated observations
        Ql = H @ np.linalg.inv(self.N) @ H.T

        # Covariance matrix of the estimated observations
        Cl = self.sigma0 ** 2 * Ql

        # Cofactor matrix of the residuals
        if not np.isfinite(np.linalg.cond(self.W)):
            log.warn("Error by computing the inverse of observation weight matrix W.")
        Qv = np.linalg.inv(W) - H @ self.Qx @ H.T

        # Redundancy
        r = np.diag(Qv @ W)

    def _print_attributes(self):

        attributes = {
            "num_obs": "number of observations",
            "num_unknowns": "number of unknowns",
            "degree_of_freedom": "degree of freedom",
            "H": "design matrix",
            "z": "observations",
            "W": "weight matrix",
            "param_names": "Parameter names",
            "x0": "apriori values of estimated parameters",
            "dx": "corrections of unknowns",
            "x_hat": "estimated solution",
            "N": "normal equation",
            "Cx": "covariance matrix of the unknowns",
            "v": "residuals",
            "sigma0": "estimated standard deviation of unit weight",
            "sigmax": "standard deviation of the unknowns",
        }

        for attribute, description in attributes.items():
            print(f"{attribute}: {description}\n{getattr(self, attribute)}\n")
