"""Get apriori data for JPL ephemerides

Description:
------------

Reads data for ephemerides using the jplephem package, which can be used to read Satellite Planet Kernel (SPK) files
provided by the Navigation and Ancillary Information Facility (NAIF) at NASA JPL. The ephemerides are read from the SPK
file 'ephemerides' listed in files.conf.

The JPL ephemeris (e.g. DE421) are aligned to the International Celestial Reference Frame (ICRF). The Barycentric
Dynamical Time (TDB) is used for accessing the JPL ephemeris (see Chapter 3 in :cite:`iers2010`).

The following convention for describing the solar objects are used by accessing SPK files, which is based on NAIF
Integer ID codes definition :cite:`naifid`:

    ==========  ============================================
     NAIF ID     Description of object
    ==========  ============================================
        0        solar system barycenter
        1        mercury barycenter
        2        venus barycenter
        3        earth barycenter   (Earth-Moon barycenter)
        4        mars barycenter
        5        jupiter barycenter
        6        saturn barycenter
        7        uranus barycenter
        8        neptune barycenter
        9        pluto barycenter
       10        sun                (Sun mass center)
      199        mercury            (Mercury mass center)
      299        venus              (Venus mass center)
      301        moon               (Moon mass center)
      399        earth              (Earth mass center)
    ==========  ============================================

The barycenter represents the center of mass of a group of bodies. For planets without moons like Mercury and Venus the
barycenter location coincides with the body center of mass. However, do not infer you may interchange use of the planet
barycenter ID and the planet ID. A barycenter has no radii, right ascension/declination of the pole axis, etc. Use the
planet ID when referring to a planet or any property of that planet.

Examples:
---------

To use the ephemerides, simply get them from the apriori-package::

    eph = apriori.get('ephemerides', time=time)

The names of objects supported by the ephemerides file can be listed::

    > eph.names
    ['earth', 'earth barycenter', 'earth moon barycenter', 'earth-moon barycenter', 'earth_barycenter', 'emb',
     'jupiter barycenter', 'jupiter_barycenter', 'mars barycenter', 'mars_barycenter', 'mercury', 'mercury barycenter',
     'mercury_barycenter', 'moon', 'neptune barycenter', 'neptune_barycenter', 'pluto barycenter', 'pluto_barycenter',
     'saturn barycenter', 'saturn_barycenter', 'solar system barycenter', 'solar_system_barycenter', 'ssb', 'sun',
     'uranus barycenter', 'uranus_barycenter', 'venus', 'venus barycenter', 'venus_barycenter']

Note that many objects have several names. For instance are `ssb`, `solar system barycenter` and
`solar_system_barycenter` different names for the same object.

To calculate a vector in either BCRS, GCRS or ITRS use the corresponding method::

    > eph.bcrs('saturn barycenter', 'venus')
    array([[  5.81353357e+11,   1.32982389e+12,   5.22346351e+11],
           [  5.80400047e+11,   1.33279492e+12,   5.23719727e+11]])

    > eph.gcrs('earth', 'earth barycenter')
    array([[-4284581.31645609, -2346398.67175642,  -714413.02172634],
           [-3684678.16574041, -3125362.44657749,  -981595.19787978]])

    > eph.itrs('earth', 'earth barycenter')
    array([[ 4882212.17477386,  -133010.97463604,  -721056.24667798],
           [ 4770255.55373369,   760435.55162305,  -987255.0146554 ]])

For convenience, a position vector of an object can be directly accessed::

    > eph.pos_bcrs('earth')     # BCRS-vector from Solar System Barycenter
    array([[ -1.26857807e+11,  -7.29571790e+10,  -3.16538026e+10],
           [ -1.25515177e+11,  -7.49544588e+10,  -3.25197196e+10]])

    > eph.pos_gcrs('moon')      # GCRS-vector from Earth Mass Center
    array([[ -3.52623481e+08,  -1.93109946e+08,  -5.87965982e+07],
           [ -3.03251110e+08,  -2.57219108e+08,  -8.07858434e+07]])

    > eph.pos_itrs('moon')      # ITRS-vector from Earth Mass Center
    array([[  4.01808840e+08,  -1.09468789e+07,  -5.93433394e+07],
           [  3.92594747e+08,   6.25842786e+07,  -8.12516495e+07]])

Similarly, velocities of objects can be computed::

    > eph.vel_bcrs('earth')
    array([[ 15327.08487282, -23240.41218388, -10075.91875001],
           [ 15752.21838668, -22991.50412836,  -9967.83854722]])

    > eph.vel_gcrs('moon')
    array([[ 484.84021811, -797.36977729, -271.44953986],
           [ 653.85202265, -680.96141269, -235.61652753]])

    > eph.vel_itrs('earth')     # Velocity is 0 because the Earth does not move
    array([[ 0.,  0.,  0.],     # relative to the Earth ...
           [ 0.,  0.,  0.]])

All vectors are in meters and meters per second (for velocity), with one vector for each time epoch specified. The
`time` argument should be a Time object so that the correct time scales can be used.

The following example shows the steps that underlies computing Earth's position referred to Solar System Barycenter::

    > from jplephem.spk import SPK
    > jd_tdb      = 2457500.5    # Julian date given in Barycentric Dynamical Time (TDB)
    > sol_bary    = 0            # NAIF ID of Solar System Barycenter
    > earth_bary  = 3            # NAIF ID of Earth-Moon Barycenter
    > earth       = 399          # NAIF ID of Earth center of mass

    # Open SPK file with planet ephemeris
    > spk = SPK.open('de430.bsp')

    # Compute Earth position in ICRF with respect to Solar System Barycenter for the given date
    > pos_earth_bary = spk[sol_bary, earth_bary].compute(jd_tdb)
    > pos_earth      = pos_earth_bary + spk[earth_bary, earth].compute(jd_tdb)
    > pos_earth
    array([ -1.26857807e+08,  -7.29571790e+07,  -3.16538026e+07])

Note that the vectors in the SPK-files are in kilometers and kilometers per day (for velocity).

TODO:
-----

Is it correct to write that JPL ephemeris are aligned to Barycentric Celestial Reference Frame (BCRS)?

If it is the case that JPL ephemeris are referred to the BCRS, also a vector between the Earth mass center and the Sun
is related to the BCRS and not to the Geocentric Celestial Reference Frame (GCRS). So far corrections between BCRS and
GCRS are not applied. What are the consequences, if we do not consider these corrections?

"""

# Standard library imports
from functools import lru_cache

# External library imports
from jplephem.spk import SPK
from jplephem import names as eph_names
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.files import dependencies
from midgard.math.unit import Unit

# Where imports
from where.lib import config


@plugins.register
def get_ephemerides(time=None, ephemerides=None):
    """Get an Ephemerides object

    The time epochs for which to calculate the ephemerides should usually be specified when getting the ephemerides. If
    they are not specified, they must be specified each time the epherides are calculated instead. The ephemerides to
    use can be specified as an input argument. If it is not specified, the config setting will be used.

    Args:
        time (Time):           Time epochs for which to calculate the ephemerides.
        ephemerides (String):  The ephemerides to use, for instance 'de430' (optional)

    Returns:
        Epherides: Object that calculates ephemerides.

    """
    ephemerides = config.tech.get("ephemerides", ephemerides).str
    return Ephemerides(time, ephemerides)


class Ephemerides:
    """A class for doing ephemerides calculations

    The jplephem package and SPK-files are used for doing the actual calculations. This class offers a more readable
    interface and caching of repeated calculations.
    """

    def __init__(self, time, ephemerides):
        """Create an Ephemerides-instance that calculates ephemerides for the given time epochs.

        It is possible to not specify the time epochs when creating the instance (set `time=None`). In this case `time`
        must be supplied as a parameter to each call calculating ephemerides.

        The SPK-file is read and parsed at the creation of this instance. In particular the names and ids of the
        available objects are read.

        Args:
            time (Time):           Time epochs for which to calculate ephemerides.
            ephemerides (String):  Name of ephemerides to use.
        """
        self.time = time
        self.ephemerides = ephemerides

        # Open the SPK-file corresponding to the ephemerides
        download_missing = config.where.files.download_missing.bool
        eph_filepath = config.files.path("ephemerides", file_vars=dict(ephemerides=ephemerides), download_missing=download_missing)
        self._spk = SPK.open(eph_filepath)  # TODO: Close file in destructor
        dependencies.add(eph_filepath, label="ephemerides")

        # Parse segments in SPK file
        self._names, self._segments = self._parse_segments()

    @property
    def names(self):
        """List names of objects available in the current ephemerides file.

        Returns:
            List:  Sorted list of strings with names of objects.
        """
        return sorted(self._names.keys())

    @lru_cache()
    def bcrs(self, from_name, to_name, time=None):
        """Calculate BCRS vectors between two objects.

        Args:
            from_name (String):  Name of object vector points from.
            to_name (String):    Name of object vector points to.
            time (Time):         Optional, time epochs for which to calculate vectors.

        Returns:
            Array: BCRS vectors, one 3-vector for each time epoch [meters].
        """
        time = self.time if time is None else time
        vector = np.zeros((3,) + time.shape)
        for segment, factor in self._find_path(from_name, to_name):
            vector += self._spk[segment].compute(time.tdb.jd) * factor

        return vector.T * Unit.kilometer2meter

    @lru_cache()
    def gcrs(self, from_name, to_name, time=None):
        """Calculate GCRS vectors between two objects.

        Currently, these are the same as the BCRS-vectors.

        Todo:
            Is it correct (enough) to use GCRS equal to BCRS or should some correction be applied?

        Args:
            from_name (String):  Name of object vector points from.
            to_name (String):    Name of object vector points to.
            time (Time):         Optional, time epochs for which to calculate vectors.

        Returns:
            Array: GCRS vectors, one 3-vector for each time epoch [meters].
        """
        return self.bcrs(from_name, to_name, time)

    @lru_cache()
    def itrs(self, from_name, to_name, time=None):
        """Calculate ITRS vectors between two objects.

        Args:
            from_name (String):  Name of object vector points from.
            to_name (String):    Name of object vector points to.
            time (Time):         Optional, time epochs for which to calculate vectors.

        Returns:
            Array: ITRS vectors, one 3-vector for each time epoch [meters].
        """
        gcrs = self.gcrs(from_name, to_name, time)

        time = self.time if time is None else time
        from where.lib import rotation

        g2i = rotation.gcrs2trs(time)

        # return position.Position(gcrs, system="gcrs", time=time).trs.val
        if time.size == 1:
            return g2i @ gcrs
        else:
            return (g2i @ gcrs[:, :, None])[:, :, 0]

    @lru_cache()
    def pos_bcrs(self, name, time=None):
        """Calculate position of object in BCRS.

        Calculated as the position of the object relative to the solar system barycenter in BCRS.

        Args:
            name (String):  Name of object.
            time (Time):    Optional, time epochs for which to calculate positions.

        Returns:
            Array: BCRS position, one 3-vector for each time epoch [meters].
        """
        return self.bcrs("solar system barycenter", name, time)

    @lru_cache()
    def pos_gcrs(self, name, time=None):
        """Calculate position of object in GCRS.

        Calculated as the position of the object relative to the earth mass center in GCRS.

        Args:
            name (String):  Name of object.
            time (Time):    Optional, time epochs for which to calculate positions.

        Returns:
            Array: GCRS position, one 3-vector for each time epoch [meters].
        """
        return self.gcrs("earth", name, time)

    @lru_cache()
    def pos_itrs(self, name, time=None):
        """Calculate position of object in ITRS.

        Calculated as the position of the object relative to the earth mass center in ITRS.

        Args:
            name (String):  Name of object.
            time (Time):    Optional, time epochs for which to calculate positions.

        Returns:
            Array: ITRS position, one 3-vector for each time epoch [meters].
        """
        return self.itrs("earth", name, time)

    @lru_cache()
    def vel(self, from_name, to_name, time=None):
        """Calculate velocity vectors between two objects.

        Args:
            from_name (String):  Name of object vector points from.
            to_name (String):    Name of object vector points to.
            time (Time):         Optional, time epochs for which to calculate velocity vectors.

        Returns:
            Array: Velocity vectors, one 3-vector for each time epoch [meters per second].
        """
        time = self.time if time is None else time
        vector = np.zeros((3,) + time.shape)
        for segment, factor in self._find_path(from_name, to_name):
            vector += self._spk[segment].compute_and_differentiate(time.tdb.jd)[1] * factor

        return vector.T * Unit.kilometer2meter / Unit.day2second

    @lru_cache()
    def vel_bcrs(self, name, time=None):
        """Calculate velocity of object in BCRS.

        Calculated as the velocity of the object relative to the solar system barycenter in BCRS.

        Args:
            name (String):  Name of object.
            time (Time):    Optional, time epochs for which to calculate velocities.

        Returns:
            Array: BCRS velocities, one 3-vector for each time epoch [meters per second].
        """
        return self.vel("solar system barycenter", name, time)

    @lru_cache()
    def vel_gcrs(self, name, time=None):
        """Calculate velocity of object in GCRS.

        Calculated as the velocity of the object relative to the earth mass center in GCRS.

        Args:
            name (String):  Name of object.
            time (Time):    Optional, time epochs for which to calculate velocities.

        Returns:
            Array: GCRS velocities, one 3-vector for each time epoch [meters per second].
        """
        return self.vel("earth", name, time)

    @lru_cache()
    def vel_itrs(self, name, time=None):
        """Calculate velocity of object in ITRS.

        Calculated as the velocity of the object relative to the earth mass center in ITRS.

        Args:
            name (String):  Name of object.
            time (Time):    Optional, time epochs for which to calculate velocities.

        Returns:
            Array: ITRS velocities, one 3-vector for each time epoch [meters per second].
        """
        vel_gcrs = self.vel_gcrs(name, time)

        time = self.time if time is None else time
        if time.isscalar:
            return time.gcrs2itrs @ vel_gcrs
        else:
            return (time.gcrs2itrs @ vel_gcrs[:, :, None])[:, :, 0]

    def _parse_segments(self):
        """Read all segments in the SPK-file

        The segments are read from the SPK-file. The names of the segments are read from the global jplephem name list,
        and gives alternative names for several objects. Any of these alternative names can be used.

        Also short aliases like `jupiter` instead of `jupiter barycenter` are added for those objects where the alias
        does not already exist. For instance (for de430), `jupiter` points to object 5, the same as `jupiter
        barycenter`, but `earth` points to object 399, while `earth barycenter` is object 3.

        These short aliases makes it easier to use for instance ephemerides and constants with the same name (like
        `eph.pos_gcrs('jupiter')` and `constant.GM_jupiter`). However, note that it is important to realize what is
        actually calculated.
        """
        segments = dict()
        for segment in self._spk.segments:
            segments[segment.center, segment.target] = 1
            segments[segment.target, segment.center] = -1

        # Add names from the global jplephem name list
        names = dict()
        for object_id in {s[0] for s in segments.keys()}:
            for name in [naif[1].lower() for naif in eph_names.target_name_pairs if naif[0] == object_id]:
                names[name] = object_id

        # Add short name aliases
        for name, object_id in names.copy().items():
            short_name = name.split()[0]
            if short_name not in names:
                names[short_name] = object_id

        return names, segments

    def _find_path(self, from_name, to_name):
        """Generate the segments necessary to calculate to find a vector between two objects

        Args:
            from_name (String):  Name of object path starts from.
            to_name (String):    Name of object path goes to.

        Returns:
            Generator: Tuples with object ids and a multiplicative factor (-1 or 1) for each segment in the path.
        """
        path = self._generate_path(self._names[from_name.lower()], self._names[to_name.lower()])
        for segment in path:
            yield segment[:: self._segments[segment]], self._segments[segment]

    def _generate_path(self, from_id, to_id, previous=None):
        """Recursively generate a path between two objects

        Effectively does a breadth first search in the `self._segments` graph.

        Args:
            from_id (Int):   Id of object path starts from.
            to_id (Int):     Id of object path goes to.
            previous_id:     Id of object previously visited (used for recursion).

        Returns:
            List: Tuples with object ids for each segment in the path.
        """
        if (from_id, to_id) in self._segments:
            return [(from_id, to_id)]
        else:
            for segment_to in [t for f, t in self._segments.keys() if f == from_id]:
                if segment_to != previous:
                    tail = self._generate_path(segment_to, to_id, previous=from_id)
                    if tail:
                        return [(from_id, segment_to)] + tail
