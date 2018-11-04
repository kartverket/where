"""A collection of Terrestrial Reference Frame sites

Description:
------------







"""

# Standard library imports
import collections

# External library imports
import numpy as np

# Where imports
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where.apriori import trf
from where.lib import util

# Position-object, TODO: replace with proper Position class
Position = collections.namedtuple("Position", ("itrs",))


class Trf(collections.UserDict):
    """A collection of Terrestrial Reference Frame sites
    """

    def __init__(self, time, reference_frames=None):
        """Set up a new Terrestrial Reference Frame object, does not parse any data

        Here we only set up the Trf-object. The individual sites are added to self.data lazily using the
        __missing__-method which is called for each key not found in self.data.

        Args:
            time (Time):                Time epochs for which to calculate positions.
            reference_frames (String):  Prioritized list of reference frames
        """
        self.time = time
        self.reference_frames = config.tech.get("reference_frames", reference_frames).list

        # Add factory for each reference frame
        self._factories = list()
        for reference_frame in self.reference_frames:
            self._factories.append(trf.get_trf_factory(time, reference_frame))

        # Cache for sites, this is populated lazily by the __missing__-method
        self.data = dict()

    def __missing__(self, key):
        """A TRF site identified by key

        """
        # Find site in the prioritized list of factories
        for factory in self._factories:
            try:
                site = factory.site(key)
                break  # Use site from this factory if it is found
            except exceptions.UnknownSiteError:
                continue  # If site is not in factory, continue to next factory
        else:  # Exit if site is not found in any factories
            log.fatal("Site '{}' not found in the reference frames {}", key, ", ".join(self.reference_frames))

        # Add site to cache
        self.data[key] = site
        return site

    @property
    def sites(self):
        """List all sites available in the reference frames
        """
        sites = set()
        for factory in self._factories:
            sites.update(factory.sites)

        return sorted(sites)

    def closest(self, pos, max_distance=5):
        """Find site closest to the given position

        Args:
            pos (Array):           3-vector with x, y and z-coordinates.
            max_distance (float):  Maximum distance around `pos` to look for a site [meter].

        Returns:
            TrfSite: Site closest to the given position. Raises `ValueError` if no site is found within `max_distance`.
        """
        distances = {k: self[k].distance_to(pos) for k in self.sites if self[k].real}
        closest = min(distances, key=distances.get)

        # Check that distance is within threshold
        if distances[closest] < max_distance:
            return self[closest]
        else:
            raise ValueError(
                "No site found within {} meters of ({:.2f}, {:.2f}, {:.2f}) in '{!r}'"
                "".format(max_distance, *pos, self)
            )

    def named_site(self, name):
        """Find site with given name
        """
        for k in self.sites:
            if name == self[k].name:
                return self[k]
        raise ValueError("No site found with name {} in '{!r}'".format(name, self))

    #
    # Dunder-methods
    #
    def __str__(self):
        return "Reference frame based on {}".format(", ".join(str(f) for f in self._factories))

    def __repr__(self):
        reference_frames = ", ".join(str(f) for f in self._factories)
        return "{}({}, '{}')".format(self.__class__.__name__, repr(self.time), reference_frames)

    def __contains__(self, item):
        if item in self.data:
            return True
        else:
            return item in self.sites

    def __iter__(self):
        for site in self.sites:
            yield self[site]


class TrfFactory:
    def __init__(self, time, version=None):
        self.name = self.__class__.__module__.split(".")[-1]
        self.time = time
        self.version = version
        self._data = None

    @property
    def data(self):
        """Data needed by this Reference Frame, lazily read by self._read_data when needed
        """
        if self._data is None:
            self._data = self._read_data()

        return self._data

    @property
    def sites(self):
        """List of all sites known by this reference frame
        """
        return self.data.keys()

    def site(self, key):
        """Positions and information about one site in the reference frame

        Args:
            key (String):  Key specifying which site to calculate position for.

        Returns:
            TrfSite:  Object with positions and information about site.
        """
        if key not in self.data:
            raise exceptions.UnknownSiteError("Unknown site '{}' in reference frame '{}'".format(key, self))

        pos_itrs = self._calculate_pos_itrs(key)
        site_info = self.data[key]
        trf_site = TrfSite(key, time=self.time, itrs=pos_itrs, source=str(self), **site_info)
        log.debug("Found site {}".format(trf_site))
        return trf_site

    #
    # Abstract methods, must be implemented by subclasses
    #
    def _read_data(self):
        """Read data needed by this Reference Frame for calculating positions of sites

        Returns:
            Dict:  Dictionary containing data about each site defined in this reference frame.
        """
        util.not_implemented()

    def _calculate_pos_itrs(self, site):
        """Calculate ITRS position of one site

        Args:
            site (String):  Key specifying which site to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 3-vector for each time epoch.
        """
        util.not_implemented()

    #
    # Dunder methods
    #
    def __str__(self):
        """String representation of factory recreates 'name:version'-string
        """
        if self.version:
            return "{}:{}".format(self.name, self.version)
        else:
            return self.name

    def __repr__(self):
        return "{}({!r}, '{}')".format(self.__class__.__name__, self.time, self.version)


class TrfSite:
    def __init__(self, key, time, itrs, source, real=True, name=None, **meta_args):
        """Constructor

        Args:
            key:     Unique identifier for the site
            time:    Epochs the site has observations for
            pos:     The position (in a TRF) of the site at each epoch
            source:  The source of the position information
            real:    Flag indicating whether the site is a real site or a dummy site
            meta:    Additional meta information
        """
        self.key = key
        self.time = time
        self.pos = Position(itrs=itrs)
        self.source = source
        self.real = real
        self.name = name
        self.meta = meta_args

        if self.name is None:
            self.name = meta_args.get(name, self.key)

    def distance_to(self, other_pos):
        """Calculate distance to some other position

        Args:
            other_pos (Array):  Other position, either as 1x3 array or nx3 with equal to length of self.time.

        Returns:
            Array:  1 or self.time.size distances.
        """
        if self.time.size == 1:
            return np.linalg.norm(self.pos.itrs - other_pos)
        return np.linalg.norm(self.pos.itrs - other_pos, axis=1).mean()

    #
    # Dunder methods
    #
    def __repr__(self):
        if self.pos.itrs.ndim == 1:
            pos = self.pos.itrs
        else:
            pos = self.pos.itrs.mean(axis=0)
        return "{}('{}', ({:.2f}, {:.2f}, {:.2f}), '{}')".format(self.__class__.__name__, self.name, *pos, self.source)
