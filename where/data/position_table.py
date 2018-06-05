"""Handling of position data inside the Where Dataset

Description:

asdf.


$Revision: 15133 $
$Date: 2018-05-18 12:50:25 +0200 (Fri, 18 May 2018) $
$LastChangedBy: dahmic $

"""

# External library imports
import numpy as np

# Where imports
from midgard.dev import console

from where.data.table import Table
from where.ext import sofa
from where.lib import cache
from where.lib import config
from where.lib import constant
from where.lib import log
from where.lib import rotation
from where.lib.unit import unit as lib_unit


class PositionTable(Table):
    """Position-klasse doc?"""
    datatype = "position"

    def __init__(self, name, num_obs, dataset):
        """Position-init doc?

            dataset:    A reference to the master Dataset. This is needed to lookup _time_obj and _other_obj.
        """
        super().__init__(name, num_obs, dataset)

        # Add units for PositionTable-properties (set by @lib_unit.register)
        self._prop_units = lib_unit.units_dict(__name__)

        # Organize data in attributes instead of a data-dict
        self._itrs = np.full((self.num_obs, 3), np.nan, dtype=float)
        self._gcrs = np.full((self.num_obs, 3), np.nan, dtype=float)
        self._time_obj = None
        self._other_obj = None
        self._dependent_on_me = list()

    def add(self, fieldname, time, itrs=None, gcrs=None, other=None, write_level=None, **_kwargs):
        """Add a field of position values

        The time must be specified when the position field is added. This can be either a TimeTable object or the
        string with the table name of a TimeTable object. In addition, one of itrs or gcrs may be specified, specifying
        the actual position vectors in either the International Terrestrial Reference System (ITRS) or the Geocentric
        Celestial Reference System (GCRS).

        Args:
            fieldnames:  String or list of strings with names of fields to be added.
            table:       String, name of table where fields are added (optional).
            time:        String, name of TimeTable with the times for the position. Needed for transformations.
            itrs:        Numpy-array with shape (num_obs, 3). ITRS position vectors.
            gcrs:        Numpy-array with shape (num_obs, 3). GCRS position vectors.
            other:       A PositionTable or DirectionTable. Used when calculating direction etc.
        """
        super().add(fieldname, write_level)
        # ToDo: Change time to time_field
        if time not in self._dataset.fields:
            raise KeyError("Time field {} does not exist in dataset".format(time))
        self._time_obj = time
        self._fields = [fieldname] + [
            fieldname + "." + a for a in self._properties if a not in ("data", "fields", "name", "num_obs")
        ]
        if other is not None:
            self.connect(other)

        if itrs is not None and gcrs is None:
            if itrs.shape != self._itrs.shape:
                raise ValueError("'itrs' must be a numpy array with shape {}".format(self._itrs.shape))
            self.itrs = itrs
        if gcrs is not None and itrs is None:
            if gcrs.shape != self._gcrs.shape:
                raise ValueError("'gcrs' must be a numpy array with shape {}".format(self._itrs.shape))
            self.gcrs = gcrs
        if itrs is not None and gcrs is not None:
            raise ValueError("At most one of 'itrs' and 'gcrs' can be set")

    def read(self, json_data, hdf5_data):
        """Read a position table from file

        The position data are read from the appropriate table in the HDF5-file and dict in the JSON-file. The time and
        other fields in the JSON-dict are simply table names that are converted to proper Table-objects later (in
        self._time_tbl and self._other_tbl) when we are certain all tables in the dataset have been read.

        Args:
            json_data:  Dict, data read from JSON-file.
            hdf5_data:  HDF5 dataset, data read from HDF5-file.

        """
        super().read(json_data, hdf5_data)
        self._itrs = hdf5_data[self.name][...]
        self._fields = json_data[self.name]["fields"]
        self._time_obj = json_data[self.name]["time"]
        self.connect(json_data[self.name]["other"])

    def write(self, json_data, hdf5_data, write_level=None):
        """Write a position table to file

        The position data dict is stored to disk by writing the itrs-positions to the HDF5-file and storing references
        to the the time- and other-objects in the JSON-file.

        Args:
            json_data:  Dict, data to be stored in JSON-file.
            hdf5_data:  HDF5 dataset, data to be stored in HDF5-file.

        """
        if not self.get_fields(write_level):
            return

        hdf5_data.create_dataset(self.name, self.itrs.shape, self.itrs.dtype)
        hdf5_data[self.name][...] = self.itrs
        json_data[self.name] = dict(
            fields=self.fields,
            time=self._time_tbl.name if self._time_tbl else None,
            other=self._other_tbl.name if self._other_tbl else None,
        )

    def copy_from(self, other_table):
        """Copy data from another position table

        Args:
            other_table:  PositionTable-object. Table to copy data from.
        """
        self._itrs = other_table.itrs.copy()
        self._fields = other_table._fields.copy()
        self._units = other_table._units.copy()
        self._write_levels = other_table._write_levels.copy()
        self._time_obj = other_table._time_tbl.name
        if other_table._other_tbl:
            self.connect(other_table._other_tbl.name)

    def subset(self, idx):
        """Remove observations from table based on idx

        Modifies the self._data numpy array in place, then updates the num_obs property.

        Args:
            idx:   Array of booleans with shape (num_obs, ). True means observation is kept, False means removed.
        """
        self._itrs = self._itrs[idx]
        self._gcrs = self._gcrs[idx]
        self._num_obs = np.sum(idx)
        self.clear_cache()

    def extend(self, other_table):
        """Add observations from another position table at the end of this table

        Args:
            other_table (PositionTable):   The other table.
        """
        self._itrs = np.vstack((self._itrs, other_table.itrs))
        self._gcrs = np.vstack((self._gcrs, np.full((other_table.num_obs, self._gcrs.shape[1]), np.nan, dtype=float)))
        self._num_obs += other_table.num_obs
        self.clear_cache()

    def plot_values(self, field):
        """Return values of a field in a form that can be plotted

        Args:
            field:   String, the field name.

        Returns:
            Numpy-array that can be plotted by for instance matplotlib.
        """
        # Plot the itrs-field by default
        if field == self.name:
            field = "{}.itrs".format(self.name)

        # Return nan's if the field is not calculable
        try:
            return getattr(self, field.split(".")[-1])
        except TypeError:
            return np.full(self.num_obs, np.nan)

    def connect(self, other):
        """Connect this PositionTable to another Position- or DirectionTable

        The other Table will be used when calculating direction, azimuth, elevation etc.

        Args:
            other:    A PositionTable or DirectionTable. Used when calculating direction etc.
        """
        if other is not None:
            self._other_obj = other if isinstance(other, str) else other.name

    def register_connection(self, other):
        self._dependent_on_me.append(other)

    def clear_cache(self, dependency="pos"):
        cache.forget_dependent_values(self, dependency)
        if dependency != "other":
            for other in self._dependent_on_me:
                other.clear_cache(dependency="other")

    def unit(self, field):
        """Unit of values in a field in table

        Args:
            field (String):  The field name.

        Returns:
            String:  The unit used for data in the field.
        """
        prop_name = field.split(".")[-1]
        return self._prop_units.get(prop_name, None)

    def __getitem__(self, key):
        """Read field data from table

        Return the appropriate list from the data dict.

        Args:
            key:   String with name of field.

        Returns:
            List of strings with field data.
        """
        if self._fields and self._fields[0].split(".")[0] == key:
            return self
        else:
            raise KeyError("{}".format(key))

    @property
    def _other_tbl(self):
        if isinstance(self._other_obj, str):
            self._other_obj = self._dataset.data[self._other_obj]
            try:
                self._other_obj.register_connection(self)
            except AttributeError:
                pass  # _other_obj does not register connections

        return self._other_obj

    @property
    def _time_tbl(self):
        if isinstance(self._time_obj, str):
            self._time_obj = self._dataset.data[self._time_obj]

        return self._time_obj

    @property
    def _time(self):
        return self._time_tbl.data

    @_time.setter
    def _time(self, value):
        self._time_obj = value

        # Force lazy recalculation of gcrs and all dependent properties
        self._gcrs[:] = None
        self.clear_cache("time")

    @property
    @lib_unit.register("meter")
    def itrs(self):
        self._assert_itrs()
        return self._itrs

    @itrs.setter
    def itrs(self, value):
        self._itrs[:] = value

        # Force lazy recalculation of gcrs and all dependent properties
        self._gcrs[:] = None
        self.clear_cache("pos")

    @property
    @lib_unit.register("meter")
    def itrs_pos(self):
        """Get position vector in ITRS

        Returns:
            numpy.ndarray:   3-dimensional array with position vectors in ITRS in [m].
        """
        return self.itrs[:, 0:3]

    def _assert_itrs(self):
        if np.any(np.isnan(self._itrs)):
            idx = np.any(np.isnan(self._itrs), axis=1)
            self._itrs[idx] = (self._time[idx].gcrs2itrs @ self._gcrs[idx, :, None])[:, :, 0]

    def add_to_itrs(self, dxyz):
        self.itrs += dxyz

    @property
    @lib_unit.register("meter")
    def gcrs(self):
        self._assert_gcrs()
        return self._gcrs

    @gcrs.setter
    def gcrs(self, value):
        self._gcrs[:] = value

        # Force lazy recalculation of itrs and all dependent properties
        self._itrs[:] = None
        self.clear_cache("pos")

    @property
    @lib_unit.register("meter")
    def gcrs_pos(self):
        """Get position vector in GCRS

        Returns:
            numpy.ndarray:   3-dimensional array with position vectors in GCRS in [m].
        """
        return self.gcrs[:, 0:3]

    def _assert_gcrs(self):
        if np.any(np.isnan(self._gcrs)):
            idx = np.any(np.isnan(self._gcrs), axis=1)
            self._gcrs[idx] = (self._time[idx].itrs2gcrs @ self._itrs[idx, :, None])[:, :, 0]

    def add_to_gcrs(self, dxyz):
        self.gcrs += dxyz

    @cache.property
    def _ref_ellipsoid(self):
        """ID of reference ellipsoid specified in the configuration file

        The IDs correspond to the ones used by the SOFA library, WGS84: 1,   GRS80: 2,   WGS72: 3.
        """
        ref_ellipsoid_cfg = config.tech.get("reference_ellipsoid")
        ref_ellipsoid = ref_ellipsoid_cfg.as_enum("reference_ellipsoid")
        log.debug("Using reference ellipsoid {} as specified in {}", ref_ellipsoid.name, ref_ellipsoid_cfg.source)
        return ref_ellipsoid

    @property
    def enu_east(self):
        """Unit vector in East direction of topocentric coordinate system

        The 1st column of the rotation matrix represents the unit vector in East direction.

        Returns:
            numpy.ndarray: Unit vector in East direction
        """
        return self._enu2itrs[:, :, 0]

    @property
    def enu_north(self):
        """Unit vector in North direction of topocentric coordinate system

        The 2nd column of the rotation matrix represents the unit vector in North direction.

        Returns:
            numpy.ndarray: Unit vector in North direction
        """
        return self._enu2itrs[:, :, 1]

    @property
    def enu_up(self):
        """Unit vector in Up direction of topocentric coordinate system

        The 3rd column of the rotation matrix represents the unit vector in Up direction.

        Returns:
            numpy.ndarray: Unit vector in Up direction
        """
        return self._enu2itrs[:, :, 2]

    def add_to_enu(self, denu):
        self.add_to_itrs(self.convert_enu_to_itrs(denu))

    @cache.dependent_property.pos
    def _enu2itrs(self):
        """Rotation matrix East, North and Up coordinate system to ITRS

        Implementation is based on Section B.2 in :cite:`subirana2013`.

        Returns:
            numpy.ndarray: Rotation matrix
        """
        lat, lon, _ = self.llh.T
        return rotation.enu2trf(lat, lon)

    @cache.dependent_property.pos
    def _itrs2enu(self):
        """Rotation matrix from ITRS to East, North and Up coordinate system

        Returns:
            numpy.ndarray: Rotation matrix
        """
        return np.transpose(self._enu2itrs, axes=[0, 2, 1])

    def convert_enu_to_itrs(self, enu):
        """Convert geocentric East, North and Up to ITRS coordinates

        The geocentric ITRS and "topocentric" coordinates (East, North, Up) are related both to the geocenter.

        Args:
            enu (numpy.ndarray):  Geocentric East, North and Up coordinates.

        Returns:
            numpy.ndarray: ITRS coordinates
        """
        return (self._enu2itrs @ enu[:, :, None])[:, :, 0]

    def convert_itrs_to_enu(self, itrs):
        """Convert ITRS to geocentric East, North and Up coordinates

        The geocentric ITRS and "topocentric" coordinates (East, North, Up) are related both to the geocenter.

        Args:
            itrs (numpy.ndarray):  ITRS coordinates

        Returns:
            numpy.ndarray: Geocentric East, North and Up coordinates
        """
        return (self._itrs2enu @ itrs[:, :, None])[:, :, 0]

    def convert_itrs_to_gcrs(self, itrs):
        return (self._time.itrs2gcrs @ itrs[:, :, None])[:, :, 0]

    def convert_gcrs_to_itrs(self, gcrs):
        return (self._time.gcrs2itrs @ gcrs[:, :, None])[:, :, 0]

    def convert_enu_to_gcrs(self, enu):
        return (self._time.itrs2gcrs @ self._enu2itrs @ enu[:, :, None])[:, :, 0]

    def convert_gcrs_to_enu(self, gcrs):
        return (self._itrs2enu @ self._time.gcrs2itrs @ gcrs[:, :, None])[:, :, 0]

    @cache.dependent_property.pos
    def llh(self):
        """Transform geocentric coordinates to geodetic using the specified reference ellipsoid

        Returns:
            numpy.ndarray: Geodetic coordinates (latitude in radians, longitude in radians, height in meters)
        """
        llh = np.empty((self.num_obs, 3))
        for obs, itrs in enumerate(self.itrs):
            longitude, latitude, height, _ = sofa.iau_gc2gd(self._ref_ellipsoid, itrs[:3])
            llh[obs, :] = (latitude, longitude, height)

        return llh

    @cache.dependent_property.pos.other.time
    @lib_unit.register("meter")
    def vector(self):
        try:
            return self._other_tbl.gcrs[:, :3] - self.gcrs[:, :3]
        except AttributeError:
            raise TypeError("No other position or direction defined. Use connect to add.") from None

    @cache.dependent_property.pos.other.time
    @lib_unit.register("meter")
    def distance(self):
        return np.linalg.norm(self.vector, axis=1)

    @cache.dependent_property.pos.other.time
    def direction(self):
        # Assume self._other_tbl is a PositionTable
        try:
            return self.vector / self.distance[:, None]
        except TypeError:
            pass

        # Assume self._other_tbl is a DirectionTable
        from where import apriori

        eph = apriori.get("ephemerides", time=self._time)
        vel = eph.vel_bcrs("earth") + self.gcrs_vel
        try:
            return (
                self._other_tbl.unit_vector
                + vel
                / constant.c
                - self._other_tbl.unit_vector
                * np.sum(self._other_tbl.unit_vector * vel, axis=1, keepdims=True)
                / constant.c
            )
        except AttributeError:
            raise TypeError("No other position or direction defined. Use connect to add.") from None

    @cache.dependent_property.pos.other.time
    @lib_unit.register("radians")
    def azimuth(self):
        east_proj = (self.direction[:, None, :] @ self._time.itrs2gcrs @ self.enu_east[:, :, None])[:, 0, 0]
        north_proj = (self.direction[:, None, :] @ self._time.itrs2gcrs @ self.enu_north[:, :, None])[:, 0, 0]
        return np.arctan2(east_proj, north_proj)

    @cache.dependent_property.pos.other.time
    @lib_unit.register("radians")
    def elevation(self):
        up_proj = (self.direction[:, None, :] @ self._time.itrs2gcrs @ self.enu_up[:, :, None])[:, 0, 0]
        return np.arcsin(up_proj)

    @cache.dependent_property.pos.other.time
    @lib_unit.register("radians")
    def zenith_distance(self):
        return np.pi / 2 - self.elevation

    def __str__(self):
        """A string describing the table

        The string describes the table, including listing the table name, datatype and all the fields. This string is
        mainly meant to be read by humans. In addition to what is shown by default by Table, we add information about
        which time and other object the position is connected to.

        Returns:
            A string describing the table.
        """
        name_and_type = "  {} ({}): ".format(self.name, self.datatype)
        fields = ", ".join(self.fields)
        time_and_other = " (time: {}, other: {})".format(
            self._time_tbl.name, self._other_tbl.name if self._other_tbl else None
        )
        return console.fill("".join((name_and_type, fields, time_and_other)), hanging=len(name_and_type))
