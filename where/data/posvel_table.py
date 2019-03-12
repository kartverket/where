"""Handling of position and velocity data inside the Where Dataset

Description:
------------

Handle both position and velocity in one structure.

-------

"""

# External library imports
import numpy as np

# Where imports
from where import apriori
from where.data.position_table import PositionTable
from where.lib import cache
from midgard.math.constant import constant
from where.lib import log
from where.lib import mathp
from where.lib import rotation
from where.lib.unit import Unit


class PosVelTable(PositionTable):
    """Function-klasse doc?"""

    datatype = "posvel"

    def __init__(self, name, num_obs, dataset):
        super().__init__(name, num_obs, dataset)

        # Add units for PosVelTable-properties (set by @Unit.register)
        self._prop_units.update(Unit.units_dict(__name__))

        # Add columns for velocity
        self._itrs = np.full((self.num_obs, 6), np.nan, dtype=float)
        self._gcrs = np.full((self.num_obs, 6), np.nan, dtype=float)

    def _assert_itrs(self):
        if np.any(np.isnan(self._itrs)):
            self._itrs[:, 0:3] = (self._time.gcrs2itrs @ self._gcrs[:, 0:3, None])[:, :, 0]
            self._itrs[:, 3:6] = (
                self._time.gcrs2itrs_dot @ self._gcrs[:, 0:3, None] + self._time.gcrs2itrs @ self._gcrs[:, 3:6, None]
            )[:, :, 0]

    def _assert_gcrs(self):
        if np.any(np.isnan(self._gcrs)):
            self._gcrs[:, 0:3] = (self._time.itrs2gcrs @ self._itrs[:, 0:3, None])[:, :, 0]
            self._gcrs[:, 3:6] = (
                self._time.itrs2gcrs_dot @ self._itrs[:, 0:3, None] + self._time.itrs2gcrs @ self._itrs[:, 3:6, None]
            )[:, :, 0]

    @property
    @Unit.register("meter per second")
    def itrs_vel(self):
        """Get velocity vector in ITRS

        Returns:
            numpy.ndarray:   3-dimensional array with velocity vectors in ITRS in [m].
        """
        return self.itrs[:, 3:6]

    @property
    @Unit.register("meter per second")
    def gcrs_vel(self):
        """Get velocity vector in GCRS

        Returns:
            numpy.ndarray:   3-dimensional array with velocity vectors in GCRS in [m].
        """
        return self.gcrs[:, 3:6]

    @cache.dependent_property.pos
    @Unit.register("radians")
    def eccentric_anomaly(self):
        r"""Compute eccentric anomaly.

        Determination of eccentric anomaly is based on Eq. 2.64 in :cite:`montenbruck2012`.

        Returns:
            numpy.ndarray: Eccentric anomaly in [rad]
        """
        E = self.kepler[:, 5]
        return E

    @cache.dependent_property.pos
    def kepler(self):
        """Compute Keplerian elements for elliptic orbit based on orbit position and velocity vector given in ITRS.

        The used equations are described in Section 2.2.4 in Montenbruck :cite:`montenbruck2012`.

        The position and velocity vector in ITRS and GM must be given in consistent units, which are [m], [m/s] and
        [m^3/s^2]. The resulting unit of the semimajor axis is implied by the unity of the position vector, i.e. [m].

        .. note::
        The function cannot be used with position/velocity vectors describing a circular or non-inclined orbit.

        Returns:
            tuple with numpy.ndarray types: Tuple with following Keplerian elements:

        ===============  ======  ==================================================================================
         Keys             Unit     Description
        ===============  ======  ==================================================================================
         a                m       Semimajor axis
         e                        Eccentricity of the orbit
         i                rad     Inclination
         Omega            rad     Right ascension of the ascending node
         omega            rad     Argument of perigee
         E                rad     Eccentric anomaly
        ===============  ======  ==================================================================================
        """
        r_norm = np.linalg.norm(self.itrs_pos, axis=1)  # Norm of position vector
        v_norm = np.linalg.norm(self.itrs_vel, axis=1)  # Norm of velocity vector

        # Normalized areal velocity (unit vector)
        h = np.cross(self.itrs_pos, self.itrs_vel)
        h_norm = np.linalg.norm(h, axis=1)
        h_unit = np.zeros((len(h), 3))
        h_unit[:, 0] = h[:, 0] / h_norm  # TODO better solution???
        h_unit[:, 1] = h[:, 1] / h_norm
        h_unit[:, 2] = h[:, 2] / h_norm

        i = np.arctan2(np.sqrt(h_unit[:, 0] ** 2 + h_unit[:, 1] ** 2), h_unit[:, 2])  # Inclination
        Omega = np.arctan2(h_unit[:, 0], -h_unit[:, 1])  # Right ascension of ascending node
        a = 1 / (2.0 / r_norm - v_norm ** 2 / constant.GM)  # Semi-major axis
        p = h_norm ** 2 / constant.GM  # Semi-latus rectum
        e = np.sqrt(1 - p / a)  # Eccentricity
        n = np.sqrt(constant.GM / a ** 3)  # Mean motion
        E = np.arctan2(
            np.einsum("ij, ij->i", self.itrs_pos, self.itrs_vel), (a ** 2 * n) * (1 - r_norm / a)  # Eccentric anomaly
        )
        vega = np.arctan2(np.sqrt(1 - e ** 2) * np.sin(E), (np.cos(E) - e))  # True anomaly
        u = np.arctan2(
            self.itrs_pos[:, 2],  # Argument of latitude
            -self.itrs_pos[:, 0] * h_unit[:, 1] + self.itrs_pos[:, 1] * h_unit[:, 0],
        )
        omega = u - vega  # Argument of perigee

        # TODO Better solution? Use of np.all() np.any()???
        for idx in range(0, len(omega)):
            if omega[idx] < 0:
                omega[idx] = omega[idx] + 2 * np.pi

        return np.vstack((a, e, i, Omega, omega, E)).T

    @cache.dependent_property.pos
    @Unit.register("radians")
    def mean_anomaly(self):
        r"""Compute mean anomaly.

        Determination of mean anomaly is based on Eq. 2.65 in :cite:`montenbruck2012`.

        Returns:
            numpy.ndarray: Mean anomaly in [rad]
        """
        kepler = self.kepler
        e = kepler[:, 1]
        E = kepler[:, 5]
        return E - e * np.sin(E)

    @cache.dependent_property.pos
    @Unit.register("radians")
    def true_anomaly(self):
        r"""Compute true anomaly.

        Determination of true anomaly is based on Eq. 2.67 in :cite:`montenbruck2012`.

        Returns:
            numpy.ndarray: True anomaly in [rad]
        """
        kepler = self.kepler
        e = kepler[:, 1]
        E = kepler[:, 5]
        return np.arctan2(np.sqrt(1 - e ** 2) * np.sin(E), (np.cos(E) - e))

    def convert_kepler_to_itrs(self, kepler):
        r"""Convert Keplerian elements to position and velocity vector in ITRS.

        Determination of mean anomaly is based on Eq. 2.65 in :cite:`montenbruck2012`.

        Args:

            kepler (numpy.ndarray): Numpy array with following Keplerian elements

        ===============  ======  ======  ==========================================================================
         Element          Index   Unit    Description
        ===============  ======  ======  ==========================================================================
         a                0       m       Semimajor axis
         e                1               Eccentricity of the orbit
         i                2       rad     Inclination
         Omega            3       rad     Right ascension of the ascending node
         omega            4       rad     Argument of perigee
         E                5       rad     Eccentric anomaly
        ===============  ======  ======  ==========================================================================

        Returns:
            numpy.ndarray: Position vector in [m] (first tree elements of the row) and velocity vector in [m/s] (last
                           three elements of the row)
        """
        itrs_pos, itrs_vel = self._keplerian_to_geocentric(
            kepler[:, 0], kepler[:, 1], kepler[:, 2], kepler[:, 3], kepler[:, 4], kepler[:, 5], anomaly="E"
        )
        return np.hstack((itrs_pos, itrs_vel))

    def _keplerian_to_geocentric(self, a, e, i, Omega, omega, x_anomaly, anomaly="E"):
        r"""Compute orbit position and velocity vector in geocentric equatorial coordinate system based on Keplerian
        elements for elliptic orbits.

        The transformation is done either for float type or numpy array/list type input arguments. The implementation
        is based on Section 2.2.3 in :cite:`montenbruck2012`.

        Args:
            a (numpy.ndarray):          Semimajor axis of orbit in [m]
            e (numpy.ndarray):          Eccentricity of the orbit
            i (numpy.ndarray):          Inclination of orbital plane in [rad]
            Omega (numpy.ndarray):      Right ascension of ascending node in [rad]
            omega (numpy.ndarray):      Argument of perigee in [rad]
            x_anomaly (numpy.ndarray):  Anomaly in [rad], which can be either mean, eccentric or true anomaly in
                                        dependency of 'anomaly' argument
            anomaly (str):              Anomaly, which can be 'M' (mean anomaly), 'E' (eccentric anomaly) or
                                        'v' (true anomaly)

        Returns:
            tuple with numpy.ndarray types:         with following elements

        ==========  =====  =======================================================================
         Element    Unit   Description
        ==========  =====  =======================================================================
         R           m     Orbit position vector in geocentric equatorial coordinate system
         V           m     Orbit velocity vector in geocentric equatorial coordinate system
        ==========  =====  =======================================================================

        """

        # Determine eccentric anomaly E
        if anomaly == "M":  # mean anomaly in [rad]
            E = get_eccentric_anomaly(x_anomaly, e)  # TODO: get_eccentric_anomaly is not defined ...
        elif anomaly == "E":  # eccentric anomaly in [rad]
            E = x_anomaly
        elif anomaly == "v":  # true anomaly in [rad]  (see Eq. 3.5 in [2])
            E = 2 * np.atan(np.tan(v / 2) * np.sqrt(1 - e / (1 + e)))  # TODO: v is not defined ...
        else:
            log.fatal(
                "Anomaly flag '{}' is wrong. It should be either 'M' (mean anomaly), 'E' (eccentric anomaly) "
                "or 'v' (true anomaly).",
                anomaly,
            )

        num_obs = len(E)
        r_orb = np.zeros((num_obs, 3))
        v_orb = np.zeros((num_obs, 3))
        R = np.zeros((num_obs, 3))
        V = np.zeros((num_obs, 3))

        cosE = np.cos(E)
        sinE = np.sin(E)
        fac = np.sqrt((1 - e) * (1 + e))
        r = a * (1 - e * cosE)  # Distance
        v = np.sqrt(constant.GM * a) / r  # Velocity

        for idx in range(0, num_obs):  # Todo: Vectorize

            # Transformation from spherical to cartesian orbital coordinate system
            r_orb = np.array([[a[idx] * (cosE[idx] - e[idx]), a[idx] * fac[idx] * sinE[idx], 0.0]])
            v_orb = np.array([[-v[idx] * sinE[idx], v[idx] * fac[idx] * cosE[idx], 0.0]])

            # Transformation from cartesian orbital to geocentric equatorial coordinate system
            PQW = rotation.R3(-Omega[idx]).dot(rotation.R1(-i[idx])).dot(rotation.R3(-omega[idx]))

            R[idx] = (PQW.dot(r_orb.T)).T
            V[idx] = (PQW.dot(v_orb.T)).T

        return R, V

    @cache.dependent_property.pos
    def _itrs2acr(self):
        """Transformation matrix from ITRS to local orbital reference system given with along-track, cross-track and
        radial.

        Returns:
            numpy.ndarray: transformation matrix with shape (num_obs, 3, 3).
        """
        r_unit = mathp.unit_vector(self.itrs_pos)  # unit vector of satellite position and in radial direction
        v_unit = mathp.unit_vector(self.itrs_vel)  # unit vector of satellite velocity
        c_unit = mathp.unit_vector(np.cross(r_unit, v_unit))  # unit vector cross-track direction
        a_unit = mathp.unit_vector(np.cross(c_unit, r_unit))  # unit vector along-track direction

        itrs2acr = np.stack((a_unit, c_unit, r_unit), axis=1)

        return itrs2acr

    @cache.dependent_property.pos
    def _acr2itrs(self):
        """Transformation matrix from local orbital reference system given with along-track, cross-track and radial
        directions to ITRS.

        Returns:
            numpy.ndarray: transformation matrix with shape (num_obs, 3, 3).
        """
        return np.transpose(self._itrs2acr, axes=[0, 2, 1])

    def convert_itrs_to_acr(self, itrs):
        """Transform ITRS coordinates into local orbital reference system (along-track, cross-track and radial)

        Returns:
            numpy.ndarray: Along-track, cross-track and radial coordinates
        """
        return (self._itrs2acr @ itrs[:, :, None])[:, :, 0]

    def convert_acr_to_itrs(self, acr):
        """Transform local orbital reference system (along-track, cross-track and radial) coordinates into ITRS
        coordinates.

        Returns:
            numpy.ndarray: ITRS coordinates.
        """
        return (self._acr2itrs @ acr[:, :, None])[:, :, 0]

    @cache.dependent_property.pos
    def itrs_pos_sun(self):
        """Determine unit vector pointing from given position to Sun in ITRS

        The determination of the vector pointing from given position to Sun is based on Eq. 5.77 in
        :cite:`subirana2013`.

        Returns:
            numpy.ndarray:  unit vectors pointing from given position to Sun in ITRS and in unit of meters
        """
        # Get Sun position vector in ITRS

        eph = apriori.get("ephemerides", time=self._time.tdb)  # TODO: is self._time.tdb correct

        # TODO:
        # Actually the JPL ephemeris are given in the BCRS with Solar System barycenter as origin and not Earth mass
        # center. So in principle the sun position vector has to be transformed from the BCRS to the GCRS. What are
        # the consequences, if we do not consider these corrections?
        earth_sun = eph.itrs("earth", "sun")  # vector pointing from Earth to Sun mass center

        # Determination of vector between given position and Sun
        pos_sun = earth_sun - self.itrs_pos
        pos_sun_unit = mathp.unit_vector(pos_sun)

        # +DEBUG: Use same model as gLAB for determination of sun-earth vector
        # glab_sun_earth = gnss.findsun(self._time.tdb)
        # glab_pos_sun = glab_sun_earth - self.itrs_pos
        # pos_sun_unit = mathp.unit_vector(glab_pos_sun)
        # -DEBUG

        return pos_sun_unit

    # TODO def gcrs_pos_sun(self):

    @cache.dependent_property.pos
    def _itrs2yaw(self):
        """Transformation matrix from ITRS to yaw-steering reference system

        The yaw-steering reference system given with x-axis lying in the Earth-Satellite-Sun plane, y-axis as the
        normal vector of the Earth-Satellite-Sun plane and the z-axis pointing to the Earth's center.

        Returns:
            numpy.ndarray: transformation matrix with shape (num_obs, 3, 3).
        """
        z_unit = -mathp.unit_vector(self.itrs_pos)  # unit vector of z-axis
        sat_sun_unit = self.itrs_pos_sun  # unit vector pointing from satellite position to Sun
        y_unit = mathp.unit_vector(np.cross(z_unit, sat_sun_unit))  # unit vector of y-axis
        x_unit = mathp.unit_vector(np.cross(y_unit, z_unit))  # unit vector of x-axis

        itrs2yaw = np.stack((x_unit, y_unit, z_unit), axis=1)

        return itrs2yaw

    @cache.dependent_property.pos
    def _yaw2itrs(self):
        """Transformation matrix from yaw-steering reference system to ITRS.

        Returns:
            numpy.ndarray: transformation matrix with shape (num_obs, 3, 3).
        """
        return np.transpose(self._itrs2yaw, axes=[0, 2, 1])

    def convert_itrs_to_yaw(self, itrs):
        """Transform ITRS coordinates into yaw-steering reference system

        Returns:
            numpy.ndarray: Yaw-steering reference system coordinates
        """
        return (self._itrs2yaw @ itrs[:, :, None])[:, :, 0]

    def convert_yaw_to_itrs(self, yaw):
        """Transform yaw-steering reference system coordinates into ITRS coordinates.

        Returns:
            numpy.ndarray: ITRS coordinates.
        """
        return (self._yaw2itrs @ yaw[:, :, None])[:, :, 0]
