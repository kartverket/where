"""Abstract class defining tables of different datatypes in a Dataset

Description:

asdf.

"""
# Standard library imports
import copy

# Where imports
from midgard.dev import console

from where.lib import cache
from where.lib import enums


class Table(object):
    datatype = None

    def __init__(self, name, num_obs, dataset):
        self._name = name
        self._num_obs = num_obs
        self._dataset = dataset
        self._data = None
        self._fields = list()
        self._units = dict()
        self._write_levels = dict()

    @classmethod
    def get_default_tablename(cls, fieldname):
        """Default name of table

        This is the table name that will be used when adding a field with this datatype unless another one is specified
        (see Dataset._add). This default implementation sets the default table name to be <fieldname> which means that
        each field will be stored in its own table. Different Tables are free to override this method, for instance to
        by default store all fields in one table with a name like <cls.datatype>_data.

        Args:
            fieldname:   String, the name of the field.

        Returns:
            String with default name of table.
        """
        return fieldname

    def add(self, fieldname, write_level=None, **_kwargs):
        """Add a field to a data table

        This method will be called once for each field added to this table from a dataset (using the
        add_datatype-method). See Dataset._add for details.

        The implementation of this method is responsible for adding the field to self._data and self._fields. It is
        also possible to add extra fields that typically will be calculated from the given field. These extra fields
        should be added to the self._fields-list as <fieldname>.<extra_field> (see PositionTable for an example of
        this).

        Note that the docstring of this function is copied to the corresponding add_datatype-method so args should
        include fieldnames (plural) and table as below, as well as any other arguments explicitly defined by the
        add-method. The docstring does not need to define a return value as that is eaten up by the
        add_datatype-method.

        Args:
            fieldnames:  String or list of strings with names of fields to be added.
            table:       String, name of table where fields are added (optional).
        """
        write_level = max(enums.get_enum("write_level")).name if write_level is None else write_level
        self._write_levels[fieldname] = enums.get_value("write_level", write_level)

    def read(self, json_data, hdf5_data):
        """Read a data table from file

        This method is called from Dataset.read.

        This method is responsible for recovering its data from the json_data and hdf5_data structures, that were
        written by write(). Just like add(), read needs to register all fields in its self._fields dictionary so that
        the fields can be calles as attributes on the Table-object.

        Args:
            json_data:  Dict, data read from JSON-file.
            hdf5_data:  HDF5 dataset, data read from HDF5-file.
        """
        self._units.update(json_data.get("_units", dict()).get(self.name, dict()))
        write_levels = json_data.get("_write_levels", dict()).get(self.name, dict())
        self._write_levels.update({f: enums.get_value("write_level", wl) for f, wl in write_levels.items()})

    def write(self, json_data, hdf5_data):
        """Write a data table to file

        This method is called from Dataset.write.

        This method is responsible for adding its data to the json_data and hdf5_data structures in such a way that the
        data can be recovered later by read(). Note that it is enough to store the data in the data structures.
        Actually storing the data to file is done by Dataset.write.

        Args:
            json_data:  Dict, data to be stored in JSON-file.
            hdf5_data:  HDF5 dataset, data to be stored in HDF5-file.
        """
        import sys

        class_name = type(self).__name__
        method_name = sys._getframe().f_code.co_name
        raise NotImplementedError(
            "'{}' class has not implemented method '{}'" "".format(class_name, method_name)
        ) from None

    def copy_from(self, other_table):
        """Copy data from another table

        This method is called from Dataset.copy_from.

        This method is reponsible for copying data from other_table to itself. It should take care to do an actual copy
        of the data, and not just copy references (i.e. it must do a proper deep copy).

        Args:
            other_table:  Table-object. Table to copy data from.
        """
        self._data = copy.deepcopy(other_table._data)
        self._fields = other_table._fields.copy()
        self._units = other_table._units.copy()
        self._write_levels = other_table._write_levels.copy()

    def subset(self, idx):
        """Remove observations from table based on idx

        This method is called from Dataset.subset.

        This method should modify the contents of self._data such that only the observations for which idx == True are
        kept. It must also update the num_obs-property (by updating self._num_obs) so that it returns the correct
        number of observations.

        Args:
            idx:   Array of booleans with shape (num_obs, ). True means observation is kept, False means removed.
        """
        import sys

        class_name = type(self).__name__
        method_name = sys._getframe().f_code.co_name
        raise NotImplementedError(
            "'{}' class has not implemented method '{}'" "".format(class_name, method_name)
        ) from None

    def extend(self, other_table):
        """Add observations from another table at the end of this table

        This method is called from Dataset.extend.

        This method should modify the contents of self._data such that it adds the contents of the `other_table` at the
        end of the table. It must also update the num_obs-property (by updating self._num_obs) so that it returns the
        correct number of observations.

        Args:
            other_table (Table):   The other table.
        """
        import sys

        class_name = type(self).__name__
        method_name = sys._getframe().f_code.co_name
        raise NotImplementedError(
            "'{}' class has not implemented method '{}'" "".format(class_name, method_name)
        ) from None

    def filter_field(self, field, filter_value):
        """Filter observations in one field

        This method is called from Dataset.filter.

        This method should return a boolean numpy-array of shape (num_obs, ). The values of the array should be True
        for each observation with value equal to filter_value, False otherwise.

        Args:
            field:        String, name of field.
            filter_value: Scalar, value to compare to.

        Returns:
            Numpy array of booleans with shape (num_obs, ).
        """
        import sys

        class_name = type(self).__name__
        method_name = sys._getframe().f_code.co_name
        raise NotImplementedError(
            "'{}' class has not implemented method '{}'" "".format(class_name, method_name)
        ) from None

    def plot_values(self, field):
        """Return values of a field in a form that can be plotted

        Args:
            field:   String, the field name.

        Returns:
            Numpy-array that can be plotted by for instance matplotlib.
        """
        return getattr(self, field.split(".")[-1])

    def as_dict(self, use_plot_values=False, fields=None):
        """Return a representation of the table as a dict

        The dictionary should also be able to be used to initialize a Pandas Dataframe.

        Args:
            use_plot_values (Boolean):  Use the plot_values instead of regular values for each field.
            fields (List):              Field names that should be included, default is to include all fields.

        Returns:
            Dict: A representation of the dataset as a dictionary.
        """
        table_dict = dict()
        fields = self.fields if fields is None else fields
        for field_name in fields:
            field = field_name.split(".")[-1]
            try:
                values = self.plot_values(field) if use_plot_values else getattr(self, field)
                if values.ndim == 1:
                    table_dict[field_name.replace(".", "_")] = values
                else:
                    for i, val in enumerate(values.T):
                        table_dict["{}_{}".format(field_name.replace(".", "_"), i)] = val
            except (AttributeError, TypeError):
                continue

        return table_dict

    @property
    def data(self):
        """Actual data stored in table

        Returns:
            Table data. The datatype will depend on the table.
        """
        return self._data

    @property
    def fields(self):
        """List names of fields in table

        Returns:
            List: Strings with names of fields.
        """
        return sorted(self._fields)

    def get_fields(self, write_level=None):
        """List names of fields in table with write level at or above `write_level`

        Args:
            write_level (String):  Lowest write level allowed for fields.

        Returns:
            List:  Strings with names of fields.
        """
        if write_level is None:
            level = min(enums.get_enum("write_level"))
        else:
            level = enums.get_value("write_level", write_level)

        return sorted(f for f, wl in self._write_levels.items() if wl >= level)

    @property
    def _write_level_strings(self):
        """List write levels of fields in table

        Returns:
            List:  Strings with write levels of fields.
        """
        return {f: wl.name for f, wl in self._write_levels.items()}

    @property
    def name(self):
        """Name of table

        Returns:
            String: name of table.
        """
        return self._name

    @property
    def num_obs(self):
        """Number of observations (= rows) in table

        Returns:
            Int: number of observations in table.
        """
        return self._num_obs

    def unit(self, field):
        """Unit of values in a field in table

        Args:
            field (String):  The field name.

        Returns:
            String:  The unit used for data in the field.
        """
        return self._units.get(field, None)

    def set_unit(self, field, unit):
        self._units[field] = unit

    @property
    def _properties(self):
        """List names of properties in table

        Dynamically looks up all properties that does not start with an underscore. Properties are defined on the
        class, not the instance, so we need to look up __class__ on self. Use issubclass to catch special properties
        like for instance those defined in lib.cache.

        One use of this property is to determine available fields at runtime as is done in PositionTable.add.

        Returns:
            List of strings with names of properties.
        """
        s_type = type(self)
        return sorted(
            k
            for k in dir(s_type)
            if issubclass(type(getattr(s_type, k)), (property, cache.property)) and not k.startswith("_")
        )

    def __getitem__(self, key):
        """Get a field from the table using dict-notation

        This method is called when a field is accessed using square brackets, e.g. table['station']. The implementation
        of __getattr__ also calls this method when accessing fields using the dot-notation.

        Args:
            key:   String, the name of the field.

        Returns:
            Field data. The datatype will depend on the table.
        """
        import sys

        class_name = type(self).__name__
        method_name = sys._getframe().f_code.co_name
        raise NotImplementedError(
            "'{}' class has not implemented method '{}'" "".format(class_name, method_name)
        ) from None

    def __getattr__(self, key):
        """Get a field from the table using attribute-notation

        This method is called when a field is accessed using a dot, e.g. table.station. This is done by forwarding the
        call to __getitem__.

        Args:
            key:   String, the name of the field.

        Returns:
            Field data. The datatype will depend on the table.
        """
        try:
            return self.__getitem__(key)
        except KeyError:
            raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, key)) from None

    def __delitem__(self, key):
        """Delete a field in the table

        This method is called by the del keyword when the field is specified using the square bracket notation, e.g.
        `del data['station']`. If the key is one of the fields we delete it from the field list.

        Args:
            key:   String, the name of the field.
        """

        if key in self._fields:
            if isinstance(self._fields, list):
                self._fields.remove(key)
            elif isinstance(self._fields, dict):
                del self._fields[key]
        else:
            raise KeyError("{} has no field '{}'".format(type(self).__name__, key))

    def __delattr__(self, key):
        """Delete a field in the table

        This method is called by the del keyword when a field is specified using the attribute notation, e.g.
        `del data.station`. If the key is one of the fields we delete it from the field list.

        Args:
            key:   String, the name of the field.
        """
        try:
            self.__delitem__(key)
        except KeyError:
            raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, key)) from None

    def __len__(self):
        """The length of the table is the number of observations

        Returns:
            Int: The number of observations in the table.
        """
        return self.num_obs

    def __dir__(self):
        """List all fields and attributes on the table

        Overridden so that fields are included in the list. This is especially useful when exploring the table
        interactively in an environment that supports tab-completion, e.g. ipython.

        Returns:
            List of strings, representing attributes and fields on the dataset.
        """
        return super().__dir__() + [f for f in self.fields if "." not in f]

    def __repr__(self):
        """A string representing the table

        The string should equal a python call that creates the table. Note that this call will only initialize the
        table, not create all fields.

        Returns:
            A string representing the table.
        """
        cls_name_full = ".".join([self.__class__.__module__.split(".")[-1], self.__class__.__name__])
        return "{cls}('{name}', {num_obs})".format(cls=cls_name_full, name=self.name, num_obs=self.num_obs)

    def __str__(self):
        """A string describing the table

        The string describes the table, including listing the table name, datatype and all the fields. This string is
        mainly meant to be read by humans.

        Returns:
            A string describing the table.
        """
        name_and_type = "  {} ({}): ".format(self.name, self.datatype)
        return console.fill(name_and_type + ", ".join(self.fields), hanging=len(name_and_type))
