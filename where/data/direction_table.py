"""Handling of direction data inside the Where Dataset

Description:
------------

asdf.

TODO: Should we use cached properties similarly to position_table?




"""

# Standard library imports

# External library imports
import numpy as np

# Where imports
from where.data.table import Table
from where.lib.unit import unit as lib_unit


class DirectionTable(Table):
    """Direction-klasse doc?"""
    datatype = "direction"

    def __init__(self, name, num_obs, dataset):
        """Direction-init doc?

            dataset:    A reference to the master Dataset. This is needed to lookup _time_obj and _other_obj.
        """
        super().__init__(name, num_obs, dataset)

        # Add units for DirectionTable-properties (set by @lib_unit.register)
        self._prop_units = lib_unit.units_dict(__name__)

    def add(self, fieldname, ra=None, dec=None, write_level=None, **_kwargs):
        """Add a field with direction values

        Args:
            fieldnames:  String or list of strings with names of fields to be added.
            table:       String, name of table where fields are added (optional).
        """
        super().add(fieldname, write_level)
        self._data = np.zeros((self.num_obs, 3))
        self._fields = [fieldname] + [
            fieldname + "." + a for a in self._properties if a not in ("data", "fields", "name", "num_obs")
        ]

        if ra is not None and dec is not None:
            self.set(right_ascension=ra, declination=dec)

        if (ra is not None and dec is None) or (ra is None and dec is not None):
            raise ValueError("Either none or both of 'ra' and 'dec' must be set")

    def read(self, json_data, hdf5_data):
        """Read a direction table from file

        The direction data are read from the appropriate dict in the JSON-file.

        Args:
            json_data:  Dict, data read from JSON-file.
            hdf5_data:  HDF5 dataset, data read from HDF5-file.
        """
        super().read(json_data, hdf5_data)
        self._data = hdf5_data[self.name][...]
        self._fields = json_data[self.name]["fields"]

    def write(self, json_data, hdf5_data, write_level=None):
        """Write a direction table to file

        The direction data dict is stored to disk by writing the itrs-directions to the HDF5-file and storing
        references to the the time- and other-objects in the JSON-file.

        Args:
            json_data:  Dict, data to be stored in JSON-file.
            hdf5_data:  HDF5 dataset, data to be stored in HDF5-file.
        """
        hdf5_data.create_dataset(self.name, self._data.shape, self._data.dtype)
        hdf5_data[self.name][...] = self._data
        json_data[self.name] = dict(fields=self._fields)

    def subset(self, idx):
        """Remove observations from table based on idx

        Modifies the self._data numpy array in place, then updates the num_obs property.

        Args:
            idx:   Array of booleans with shape (num_obs, ). True means observation is kept, False means removed.
        """
        self._data = self._data[idx]
        self._num_obs = np.sum(idx)

    def extend(self, other_table):
        """Add observations from another direction table at the end of this table

        Args:
            other_table (DirectionTable):   The other table.
        """
        self._data = np.vstack((self._data, other_table.data))
        self._num_obs += other_table.num_obs

    def set(self, right_ascension, declination):
        self._data[:] = np.array(
            (
                np.cos(declination) * np.cos(right_ascension),
                np.cos(declination) * np.sin(right_ascension),
                np.sin(declination),
            )
        ).T

    def unit(self, field):
        """Unit of values in a field in table

        Args:
            field (String):  The field name.

        Returns:
            String:  The unit used for data in the field.
        """
        prop_name = field.split(".")[-1]
        return self._prop_units.get(prop_name, None)

    @property
    def unit_vector(self):
        norms = np.linalg.norm(self._data, axis=1, keepdims=True)
        norms[norms == 0] = 1
        return self._data / norms

    @property
    @lib_unit.register("radians")
    def declination(self):
        return np.arcsin(self.unit_vector[:, 2])

    @property
    @lib_unit.register("radians")
    def right_ascension(self):
        return np.arctan2(self.unit_vector[:, 1], self.unit_vector[:, 0])

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
