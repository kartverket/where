"""A collection of Celestial Reference Frame sources

Description:
------------

"""

# Standard library imports
import collections

# Where imports
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where.apriori import crf
from where.lib import util

# Position-object, TODO: replace with proper Position class
# First element right ascension, second element declination
# Position = collections.namedtuple("Position", ("crs",))


class Crf(collections.UserDict):
    """A collection of Celestial Reference Frame sources
    """

    def __init__(self, time, celestial_reference_frames=None):
        """Set up a new Celestial Reference Frame object, does not parse any data

        Here we only set up the Crf-object. The individual sources are added to self.data lazily using the
        __missing__-method which is called for each key not found in self.data.

        Args:
            celestial_reference_frames (String):  Prioritized list of celestial reference frames
        """
        super().__init__()
        self.time = time
        self.celestial_reference_frames = config.tech.get(
            "celestial_reference_frames", celestial_reference_frames
        ).list

        # Add factory for each celestial reference frame
        self._factories = list()
        for celestial_reference_frame in self.celestial_reference_frames:
            self._factories.append(crf.get_crf_factory(time, celestial_reference_frame))

    def __missing__(self, key):
        """A CRF source identified by key

        """
        # Find source in the prioritized list of factories
        for factory in self._factories:
            try:
                source = factory.source(key)
                break  # Use source from this factory if it is found
            except exceptions.UnknownRadioSourceError:
                continue  # If source is not in factory, continue to next factory
        else:  # Exit if site is not found in any factories
            log.fatal(
                f"Source {key!r} not found in the celestial reference frames "
                f"{', '.join(self.celestial_reference_frames)}"
            )

        # Add source to cache
        self.data[key] = source
        return source

    @property
    def sources(self):
        """List all sites available in the reference frames
        """
        sources = set()
        for factory in self._factories:
            sources.update(factory.sources)

        return sorted(sources)

    #
    # Dunder-methods
    #
    def __str__(self):
        return "Celestial reference frame based on {}".format(", ".join(str(f) for f in self._factories))

    def __repr__(self):
        celestial_reference_frames = ", ".join(str(f) for f in self._factories)
        return "{}('{}')".format(self.__class__.__name__, celestial_reference_frames)

    def __contains__(self, item):
        if item in self.data:
            return True
        else:
            return item in self.sources

    def __iter__(self):
        for source in self.sources:
            yield self[source]


class CrfFactory:
    def __init__(self, time):
        self.name = self.__class__.__module__.split(".")[-1]
        self.time = time
        self._data = None

    @property
    def data(self):
        """Data needed by this Celestial Reference Frame, lazily read by self._read_data when needed
        """
        if self._data is None:
            self._data = self._read_data()

        return self._data

    @property
    def sources(self):
        """List of all sources known by this celestial reference frame
        """
        return self.data.keys()

    def source(self, key):
        """Positions and information about one source in the celestial reference frame

        Args:
            key (String):  Key specifying which source to calculate position for.

        Returns:
            CrfSource:  Object with positions and information about source.
        """
        if key not in self.data:
            raise exceptions.UnknownRadioSourceError(
                "Unknown source '{}' in celestial reference frame '{}'".format(key, self)
            )

        pos_crs = self._calculate_pos_crs(key)
        source_info = self.data[key]

        crf_source = CrfSource(key, pos_crs, source=str(self), **source_info)
        log.debug(f"Found source {crf_source}")
        return crf_source

    #
    # Abstract methods, must be implemented by subclasses
    #
    def _read_data(self):
        """Read data needed by this Celestial Reference Frame for calculating positions of sources

        Returns:
            Dict:  Dictionary containing data about each source defined in this celestial reference frame.
        """
        util.not_implemented()

    def _calculate_pos_crs(self, source):
        """Calculate CRS position of one source

        Args:
            source (String):  Key specifying which source to calculate position for, must be key in self.data.

        Returns:
            Array:  Positions, one 2-vector for each time epoch.
        """
        util.not_implemented()

    #
    # Dunder methods
    #
    def __str__(self):
        """String representation of factory recreates 'crf1'-string
        """
        return self.name

    def __repr__(self):
        return "{}()".format(self.__class__.__name__)


class CrfSource:
    def __init__(self, key, direction, source, **meta_args):
        self.key = key
        self.pos = direction
        self.source = source
        self.meta = meta_args

        # Try to set name from meta, with fallback to key
        self.name = self.meta.get("name", self.key)

    #
    # Dunder methods
    #
    def __repr__(self):
        return "{}('{}', ({:.2f}, {:.2f}), '{}')" "".format(
            self.__class__.__name__, self.name, *self.pos.gcrs, self.source
        )
