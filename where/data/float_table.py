"""Handling of float data inside the Where Dataset

Description:

asdf.



"""

# External library imports
import numpy as np

# Where imports
from where.data.table import Table


class FloatTable(Table):
    """Float-klasse doc?"""
    datatype = "float"

    def __init__(self, name, num_obs, dataset):
        """Float-init doc?"""
        super().__init__(name, num_obs, dataset)
        self._data = np.zeros((num_obs, 0))
        self._fields = dict()

    @classmethod
    def get_default_tablename(cls, fieldname):
        """Default name of table

        Float data are by default stored in one table (numpy-array) called 'float_data'.

        Returns:
            String with default name of table.
        """
        return "{}_{}".format(cls.datatype, "data")

    def add(self, fieldname, shape=(), val=None, unit=None, write_level=None, **_kwargs):
        """Add a field of float values

        Args:
            fieldnames:  String or list of strings with names of fields to be added.
            table:       String, name of table where fields are added (optional).
            shape:       Tuple of ints specifying the shape of each observation.
        """
        super().add(fieldname, write_level)
        size = np.prod(shape, dtype=int)
        empty_field = np.zeros((self.num_obs, size))
        idx_first = self._data.shape[1]
        self._fields[fieldname] = slice(idx_first, idx_first + size), shape
        self._data = np.append(self._data, empty_field, axis=1)

        if val is not None:
            val = np.array(val)
            if val.shape != self[fieldname].shape:
                raise ValueError("The shape of 'val' must be {}, not {}".format(self[fieldname].shape, val.shape))
            self[fieldname][:] = val

        if unit is not None:
            self._units[fieldname] = unit

    def read(self, json_data, hdf5_data):
        super().read(json_data, hdf5_data)
        self._data = hdf5_data[self.name][...]
        fields = hdf5_data[self.name].attrs["fields"].split(",")
        idxes = [tuple(int(i) for i in idx.split("-")) for idx in hdf5_data[self.name].attrs["idxes"].split(",")]
        shapes = [
            tuple(int(s) for s in (shp.split("-") if shp else ()))
            for shp in hdf5_data[self.name].attrs["shapes"].split(",")
        ]
        self._fields = {f: (slice(*i), s) for f, i, s in zip(fields, idxes, shapes)}

        return fields

    def write(self, json_data, hdf5_data, write_level=None):
        fields = self.get_fields(write_level)

        # Reshape the data if not all fields are stored
        if set(self.fields) - set(fields):
            field_data = list()
            field_spec = dict()
            total_size = 0
            for field in fields:
                idx, shape = self._fields[field]
                size = np.prod(shape, dtype=int)
                field_spec[field] = slice(total_size, total_size + size), shape
                field_data.append(self._data[:, idx])
                total_size += size
            data = np.concatenate(field_data, axis=1)
        else:
            data = self._data
            field_spec = self._fields

        # Store data and some attributes/metadata in HDF5
        hdf5_data.create_dataset(self.name, data.shape, data.dtype)
        hdf5_data[self.name][...] = data

        idxes = ["{s.start}-{s.stop}".format(s=field_spec[f][0]) for f in fields]
        shapes = ["-".join(str(s) for s in field_spec[f][1]) for f in fields]
        hdf5_data[self.name].attrs["fields"] = ",".join(fields)
        hdf5_data[self.name].attrs["idxes"] = ",".join(idxes)
        hdf5_data[self.name].attrs["shapes"] = ",".join(shapes)

    def subset(self, idx):
        """Remove observations from table based on idx

        Modifies the self._data numpy array in place, then updates the num_obs property.

        Args:
            idx:   Array of booleans with shape (num_obs, ). True means observation is kept, False means removed.
        """
        self._data = self._data[idx]
        self._num_obs = self._data.shape[0]

    def extend(self, other_table):
        """Add observations from another float table at the end of this table

        Fields in only one of the tables are handled by filling with zeros.

        Args:
            other_table (FloatTable):   The other table.
        """
        # Add fields that are in other table, but not in this one
        new_fields = set(other_table.fields) - set(self.fields)
        for new_field in new_fields:
            self.add(new_field, shape=other_table[new_field].shape[1:], unit=other_table.unit(new_field))

        # Extend all fields with the needed number of observations
        empty_table = np.zeros((other_table.num_obs,) + self._data.shape[1:])
        self._data = np.vstack((self._data, empty_table))
        self._num_obs += other_table.num_obs

        # Copy data from the other table (one field at a time in case indices are different)
        for field in other_table.fields:
            self[field][-other_table.num_obs:] = other_table[field]

    def __getitem__(self, key):
        """Read field data from table

        Look up location and shape of field data in self._fields. Then return data after they are reshaped.

        Args:
            key:   String with name of field.

        Returns:
            Numpy-array with field data.
        """
        idx, shape = self._fields[key]
        return self._data[:, idx].reshape(-1, *shape)
