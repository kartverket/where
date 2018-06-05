"""Handling of time data inside the Where Dataset

Description:

asdf.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""

# External library imports
import numpy as np

# Where imports
from where.data.table import Table
from where.lib import time


class TimeTable(Table):
    """Time-klasse doc?"""
    datatype = "time"

    def add(self, fieldname, val, scale, write_level=None, **kwargs):
        super().add(fieldname, write_level)
        self._data = time.Time(val, scale=scale, **kwargs)
        self._fields = [fieldname]
        self._units[fieldname] = scale  # Use time scale as unit, mainly for visualizing time scale

    def read(self, json_data, hdf5_data):
        super().read(json_data, hdf5_data)
        self._fields = json_data[self.name]["fields"]
        scale = json_data[self.name].get("scale", "tai")  # TODO: 'tai' used earlier, `get` for backwards compatibility
        #       should be `json_data[self.name]['scale']` but will
        #       cause old datasets to become unreadable
        time_values = hdf5_data[self.name][...]
        self._data = time.Time(val=time_values[:, 0], val2=time_values[:, 1], format="jd", scale=scale)

    def write(self, json_data, hdf5_data, write_level=None):
        if not self.get_fields(write_level):
            return

        hdf5_data.create_dataset(self.name, (self.num_obs, 2), float)
        hdf5_data[self.name][:, 0] = self._data.jd1
        hdf5_data[self.name][:, 1] = self._data.jd2

        json_data[self.name] = dict(fields=self._fields, scale=self._data.scale)

    def subset(self, idx):
        self._data = self._data[idx]
        self._num_obs = np.sum(idx)

    def extend(self, other_table):
        """Add observations from another time table at the end of this table

        Fields in only one of the tables are handled by filling with empty strings.

        TODO: Make sure time scale is handled correctly

        Args:
            other_table (TimeTable):   The other table.
        """
        self._data = time.Time(
            np.hstack((self._data.utc.jd1, other_table.data.utc.jd1)),
            np.hstack((self._data.utc.jd2, other_table.data.utc.jd2)),
            format="jd",
            scale="utc",
        )
        self._num_obs += other_table.num_obs

    def filter_field(self, field, filter_value):
        """Filter observations in one field

        Args:
            field:         String, name of field.
            filter_value:  Time, value to compare to.

        Returns:
            Numpy array of booleans with shape (num_obs, ).
        """
        return np.isclose(self._data.jd, filter_value.jd, rtol=0)

    def plot_values(self, field):
        """Return values of a field in a form that can be plotted

        Args:
            field:   String, the field name.

        Returns:
            Numpy-array that can be plotted by for instance matplotlib.
        """
        return super().plot_values(field).datetime

    def as_dict(self, use_plot_values=False, fields=None):
        """Return a representation of the table as a dict

        The dictionary should also be able to be used to initialize a Pandas Dataframe.

        Args:
            use_plot_values (Boolean):  Use the plot_values instead of regular values for each field.
            fields (List):              Field names that should be included, default is to include all fields.

        Returns:
            Dict: A representation of the table as a dictionary.
        """
        if fields is None or self.name in fields:
            return {self.name: self.plot_values(self.name) if use_plot_values else self._data.datetime}
        else:
            return dict()

    def __getitem__(self, key):
        """Read field data from table

        Look up location and shape of field data in self._fields. Then return data after they are reshaped.

        Args:
            key:   String with name of field.

        Returns:
            Numpy-array with field data.
        """
        if key == self.name:
            return self._data

        # We should not end up here ...
        from where.lib import log

        log.error("Somehow key is {} but my name is {}", key, self.name)

    def __dir__(self):
        """List all fields and attributes on the table

        Overridden so that the attributes of the Time object itself is included in the list.

        Returns:
            List of strings, representing attributes and fields on the dataset.
        """
        return super().__dir__() + dir(self._data)
