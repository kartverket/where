"""Handling of text data inside the Where Dataset

Description:

asdf.



"""

# Standard library imports
from collections import UserList

# External library imports
import numpy as np

# Where imports
from where.data.table import Table


class TextList(UserList):

    def __getitem__(self, key):
        """Add fancy boolean indexing to the list
        """
        if isinstance(key, np.ndarray):
            return self.__class__([sat for idx, sat in zip(key, self.data) if idx])
        else:
            return super().__getitem__(key)

    def __setitem__(self, key, value):
        """Add fancy boolean indexing to the list
        """
        if isinstance(key, np.ndarray):
            for idx in np.where(key)[0]:
                self.data[idx] = value
        else:
            super().__setitem__(key, value)


class TextTable(Table):
    """Text-klasse doc?"""
    datatype = "text"

    def __init__(self, name, num_obs, dataset):
        """Text-init doc?"""
        super().__init__(name, num_obs, dataset)
        self._data = dict()

    @classmethod
    def get_default_tablename(cls, fieldname):
        """Default name of table

        Text data are by default stored in one table (dict) called 'text_data'.

        Returns:
            String with default name of table.
        """
        return "{}_{}".format(cls.datatype, "data")

    def add(self, fieldname, val, unit="", write_level=None, **_kwargs):
        """Add a field of text values

        Args:
            fieldnames:  String or list of strings with names of fields to be added.
            table:       String, name of table where fields are added (optional).
        """
        super().add(fieldname, write_level)
        self._data[fieldname] = np.array(val, dtype=str)
        self._fields.append(fieldname)
        self._units[fieldname] = unit

        if val is not None:
            if len(val) != self.num_obs:
                raise ValueError("'val' must be a list with length {}".format(self.num_obs))

    def read(self, json_data, hdf5_data):
        """Read a text table from file

        The text data are read from the appropriate dict in the JSON-file.

        Args:
            json_data:  Dict, data read from JSON-file.
            hdf5_data:  HDF5 dataset, data read from HDF5-file.

        Returns:
            List of strings with all field names read from file.
        """
        super().read(json_data, hdf5_data)
        self._data = {f: np.array(d, dtype=str) for f, d in json_data[self.name].items()}
        self._fields = sorted(self._data)

    def write(self, json_data, hdf5_data, write_level=None):
        """Write a text table to file

        The text data dict is simply stored as a dict in the JSON-file.

        TODO: Store text data in HDF5-file

        Args:
            json_data:  Dict, data to be stored in JSON-file.
            hdf5_data:  HDF5 dataset, data to be stored in HDF5-file.
        """
        json_data[self.name] = {f: self._data[f].tolist() for f in self.get_fields(write_level)}

    def copy_from(self, other_table):
        """Copy data from another text table

        Args:
            other_table:  TextTable-object. Table to copy data from.
        """
        super().copy_from(other_table)
        self._data = {f: v.copy() for f, v in other_table.data.items()}

    def subset(self, idx):
        """Remove observations from table based on idx

        Modifies the self._data numpy array in place, then updates the num_obs property.

        Args:
            idx:   Array of booleans with shape (num_obs, ). True means observation is kept, False means removed.
        """
        for key, fielddata in self._data.items():
            self._data[key] = fielddata[idx]
        self._num_obs = np.sum(idx)

    def extend(self, other_table):
        """Add observations from another text table at the end of this table

        Fields in only one of the tables are handled by filling with empty strings.

        Args:
            other_table (TextTable):   The other table.
        """
        fields = set(self.fields) | set(other_table.fields)
        for field in fields:
            other_data = other_table[field] if field in other_table.fields else [""] * other_table.num_obs
            if field not in self.fields:
                self._dataset.add_text(field, val=other_data)
            else:
                self._data[field] = np.concatenate((self._data[field], other_data))

        self._num_obs += other_table.num_obs

    def filter_field(self, field, filter_value):
        """Filter observations in one field

        Args:
            field:         String, name of field.
            filter_value:  String, value to compare to.

        Returns:
            Numpy array of booleans with shape (num_obs, ).
        """
        return self._data[field] == filter_value

    def plot_values(self, field):
        """Return values of a field in a form that can be plotted

        Args:
            field:   String, the field name.

        Returns:
            Numpy-array that can be plotted by for instance matplotlib.
        """
        values = super().plot_values(field)
        unique_values = sorted(set(values))
        return np.array([unique_values.index(v) + 1 for v in values])

    def as_dict(self, use_plot_values=False, fields=None):
        """Return a representation of the table as a dict

        The dictionary should also be able to be used to initialize a Pandas Dataframe.

        Args:
            use_plot_values (Boolean):  Use the plot_values instead of regular values for each field.
            fields (List):              Field names that should be included, default is to include all fields.

        Returns:
            Dict: A representation of the dataset as a dictionary.
        """
        fields = self.fields if fields is None else fields
        if use_plot_values:
            return {f: self.plot_values(f) for f in self.fields if f in fields}
        else:
            return {f: val for f, val in self._data.items() if f in fields}

    def __getitem__(self, key):
        """Read field data from table

        Return the appropriate list from the data dict.

        Args:
            key:   String with name of field.

        Returns:
            List of strings with field data.
        """
        return self._data[key]
