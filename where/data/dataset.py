"""A Dataset that handles all data in Where

Description:
------------

asdf.




"""

# Standard library imports
import copy
from datetime import date
import json

# External library imports
import numpy as np
import pandas as pd

# Where imports
from midgard.dev import console

import where
from where.data import _data
from where.lib import config
from where.lib.exceptions import FieldExistsError, InitializationError
from where.lib import files
from where.lib import log
from where.lib.time import Time


class Dataset(object):
    r"""A Dataset represents data for one given model run.

    Data are represented as fields, which are stored in tables under the hood. Each table can handle one specific data
    type. See the Table abstract class and its subclass implementations for details about each table type.

    TODO:
           - Documentation for fields (descriptions)
           - Copy from dataset
           - Dataset-id
           - Documentation

    Here is a simple sphinx math test, :math:`\frac12 \int_0^1 a^2 + b^2 \, \mathrm da`, how about that?

    And this is neat:

    .. math::
       \mathrm{correction} = \frac{- \hat K \cdot \vec b \bigl[ 1 - \frac{(1 + \gamma) U}{c^2}
                             - \frac{| \vec V_\oplus |^2}{2 c^2} - \frac{\vec V_\oplus \cdot \vec w_2}{c^2} \bigr]
                             - \frac{\vec V_\oplus \cdot \vec b}{c} \bigl[ 1
                             + \frac{\hat K \cdot \vec V_\oplus}{2 c} \bigr]}{1
                             + \frac{\hat K \cdot (\vec V_\oplus + \vec w_2)}{c}}

    Regular text?
    """

    def __init__(self, rundate, tech, stage, dataset_name, dataset_id, empty=False, **kwargs):
        """Create a new Dataset or read an existing one

        Note:
            Be aware that the implementation is dependent on ``self._fields`` being the first attribute to be set. See
            :func:`__setattr__` for more information.

        Args:
            rundate:      Date, the model run date.
            tech:         String, the technique.
            stage:        String, the stage.
            dataset_name: String, the name of the dataset.
            dataset_id:   Int, id of the dataset.
            empty:        Boolean, if False (default) will read dataset from disk if available.

        """
        self._fields = dict()
        self._data = dict()
        self._num_obs = 0
        self._default_field_suffix = None
        self._kwargs = kwargs
        self._kwargs.setdefault("session", dataset_name)  # TODO: Can this be removed?
        self.vars = dict(
            config.program_vars(
                **dict(
                    kwargs,
                    rundate=rundate,
                    tech_name=tech,
                    stage=stage,
                    dataset_name=dataset_name,
                    dataset_id=str(dataset_id),
                )
            )
        )
        self.vars.update(**kwargs)
        self.rundate = rundate
        dataset_id = _data.parse_dataset_id(rundate, tech, stage, dataset_name, dataset_id)
        self.name = "{name}/{id:04d}".format(name=dataset_name, id=dataset_id)
        self.meta = dict()

        # Try to read dataset from disk unless explicitly told to create an empty dataset
        if not empty:
            try:
                self.read()
            except FileNotFoundError:
                pass

    @classmethod
    def anonymous(cls, num_obs=None):
        """Create an anonymous dataset

        The dataset has no meaningful values for rundate, tech, stage, dataset_name, and dataset_id. These should be
        set with `rename` or `write_as` before/as the dataset is written to disk.

        Args:
            num_obs (Int):  Number of observations in dataset (optional).
        """
        dset = cls(rundate=date.today(), tech="", stage="", dataset_name="", dataset_id=0, empty=True)
        if num_obs is not None:
            dset.num_obs = num_obs

        return dset

    def read(self):
        """Read a dataset from file

        A dataset is stored on disk in two files, one JSON-file and one HDF5-file. Typically the HDF5-file is great for
        handling numeric data, while JSON is more flexible. The actual reading of the data is handled by the individual
        datatype table-classes. The dispatch to the correct class is done by functions defined in the
        :func:`Dataset._read`-method which is called by :mod:`where.data._data` when Dataset is first imported.
        """
        # Open and read JSON-file
        json_path = files.path("dataset_json", file_vars=self.vars)
        with files.open_path(json_path, mode="rt", write_log=False) as f_json:
            json_all = json.load(f_json)
        if self.name not in json_all:
            raise FileNotFoundError("Dataset {} not found in file {}".format(self.name, json_path))

        log.debug(f"Read dataset {self.vars['tech']}-{self.vars['stage']} from disk at {json_path.parent}")
        json_data = json_all[self.name]
        self._num_obs = json_data["_num_obs"]
        tables = json_data["_tables"]

        # Open HDF5-file
        with files.open_datafile("dataset_hdf5", file_vars=self.vars, mode="r", write_log=False) as f_hdf5:
            hdf5_data = f_hdf5[self.name]

            # Read data for each table by dispatching to read function based on datatype
            for table, dtype in tables.items():
                read_func = getattr(self, "_read_" + dtype)
                read_func(table, json_data, hdf5_data)

        # Add meta and vars properties
        self.meta = json_data.get("_meta", dict())
        self.vars = json_data.get("_vars", self.vars)

    def rename(self, rundate=None, tech=None, stage=None, dataset_name=None, dataset_id=None, **kwargs):
        """Rename a dataset

        Renames the dataset. In particular, this means that if the dataset is written to file it will be written to a
        different file (or different place in the same file). All arguments are optional. If they are not given, they
        keep their existing value.

        Args:
            rundate:      Date, the model run date.
            tech:         String, the technique.
            stage:        String, the stage.
            dataset_name: String, the name of the dataset.
            dataset_id:   Int, id of the dataset.
        """
        # Set rundate
        if rundate is not None:
            self.rundate = rundate

        # Use existing names as default
        tech = self.vars["tech"] if tech is None else tech
        stage = self.vars["stage"] if stage is None else stage
        dataset_name = self.dataset_name if dataset_name is None else dataset_name
        if dataset_id is None:
            dataset_id = self.dataset_id
        else:
            dataset_id = _data.parse_dataset_id(self.rundate, tech, stage, dataset_name, dataset_id)

        # Update names
        self.name = "{name}/{id:04d}".format(name=dataset_name, id=dataset_id)
        kwargs.setdefault("session", dataset_name)
        self.vars.update(dict(tech=tech, stage=stage, dataset_name=dataset_name, dataset_id=dataset_id, **kwargs))

    def write(self, write_level=None):
        """Write a dataset to file

        A dataset is stored on disk in two files, one JSON-file and one HDF5-file. Typically the HDF5-file is great for
        handling numeric data, while JSON is more flexible. The actual writing of the data is handled by the individual
        datatype table-classes. These classes are free to choose how they divide the data between the JSON- and
        HDF5-files, as long as they are able to recover all the data.
        """
        json_path = files.path("dataset_json", file_vars=self.vars)
        log.debug(f"Write dataset {self.vars['tech']}-{self.vars['stage']} to disk at {json_path.parent}")

        # Read write level from config
        write_level = config.tech.get("write_level", value=write_level).as_enum("write_level").name

        # Read existing data in JSON-file
        try:
            with files.open_path(json_path, mode="rt", write_log=False) as f_json:
                json_all = json.load(f_json)
        except FileNotFoundError:
            json_all = dict()
        json_all.setdefault(self.name, dict())
        json_data = json_all[self.name]

        # Figure out which tables have data
        tables = [t for t in self._data.values() if t.get_fields(write_level)]

        # Open HDF5-file
        with files.open_datafile("dataset_hdf5", file_vars=self.vars, mode="a", write_log=False) as f_hdf5:
            if self.name in f_hdf5:
                del f_hdf5[self.name]
            hdf5_data = f_hdf5.create_group(self.name)

            # Write data for each table (HDF5-data are automatically written to disk)
            for table in tables:
                table.write(json_data, hdf5_data, write_level)

        # Store metadata in JSON-data
        json_data["_version"] = where.__version__
        json_data["_num_obs"] = self.num_obs
        json_data["_tables"] = {tbl.name: tbl.datatype for tbl in tables}
        json_data["_units"] = {tbl.name: tbl._units for tbl in tables}
        json_data["_write_levels"] = {tbl.name: tbl._write_level_strings for tbl in tables}
        json_data["_meta"] = self.meta
        json_data["_vars"] = self.vars

        # Store last dataset_id written to
        json_all.setdefault(self.dataset_name, dict())["_last_dataset_id"] = self.dataset_id

        # Write JSON-data to file
        with files.open_path(json_path, mode="wt", write_log=False) as f_json:
            json.dump(json_all, f_json)

    def write_as(
        self, rundate=None, tech=None, stage=None, dataset_name=None, dataset_id=None, write_level=None, **kwargs
    ):
        """Write a dataset to file with a new name

        Renames the dataset and writes it to file.

        Args:
            rundate:      Date, the model run date.
            tech:         String, the technique.
            stage:        String, the stage.
            dataset_name: String, the name of the dataset.
            dataset_id:   Int, id of the dataset.
        """
        self.rename(rundate, tech, stage, dataset_name, dataset_id, **kwargs)
        self.write(write_level=write_level)

    def delete_from_file(self, tech=None, stage=None, dataset_name=None, dataset_id=None):
        """Delete this or related datasets from file

        Specify arguments relative to this dataset to find datasets which will be deleted.
        """
        # Use existing names as default
        tech = self.vars["tech"] if tech is None else tech
        stage = self.vars["stage"] if stage is None else stage
        dataset_name = self.dataset_name if dataset_name is None else dataset_name
        if dataset_id is None:
            dataset_id = self.dataset_id
        else:
            dataset_id = _data.parse_dataset_id(self.rundate, tech, stage, dataset_name, dataset_id)

        dataset_id = {dataset_id} if isinstance(dataset_id, (float, int)) else set(dataset_id)
        ids_to_delete = dataset_id & set(_data.list_dataset_ids(self.rundate, tech, dataset_name, stage, dataset_name))
        if not ids_to_delete:
            return

        # Open JSON and HDF5 file and remove datasets
        file_vars = dict(self.vars, tech=tech, stage=stage)
        json_path = files.path("dataset_json", file_vars=file_vars)
        with files.open_path(json_path, mode="rt", write_log=False) as f_json:
            json_all = json.load(f_json)
        with files.open_datafile("dataset_hdf5", file_vars=file_vars, mode="a", write_log=False) as f_hdf5:
            for id_to_delete in ids_to_delete:
                name = "{name}/{id:04d}".format(name=dataset_name, id=id_to_delete)
                del json_all[name]
                del f_hdf5[name]
                log.debug(f"Deleted {name} from dataset {tech}-{stage} at {json_path.parent}"),

        with files.open_path(json_path, mode="wt", write_log=False) as f_json:
            json.dump(json_all, f_json)

        # Delete files if all datasets are deleted
        if not any(["/" in k for k in json_all.keys()]):
            json_path.unlink()
            files.path("dataset_hdf5", file_vars=file_vars).unlink()

    def copy_from(self, other_dataset, meta_key=None):
        """Copy observations from another dataset

        Args:
            other_dset (Dataset):  The other dataset.
            meta_key (str):        Dictionary key for introduction of an additional level in dictionary.
        """
        # Check and update number of observations
        if self.num_obs and self.num_obs != other_dataset.num_obs:
            raise ValueError(
                "The other dataset '{}' has {} observations, which is incompatible with self.num_obs = {}"
                "".format(other_dataset.name, other_dataset.num_obs, self.num_obs)
            )
        self._num_obs = other_dataset.num_obs

        # Clear any existing fields and tables
        self._fields = dict()
        self._data = dict()

        # Update meta information
        if meta_key is None:
            self.meta = copy.deepcopy(other_dataset.meta)
        else:
            self.meta.setdefault(meta_key, dict())
            self.meta[meta_key] = copy.deepcopy(other_dataset.meta)

        # Create each table, copy the data and register fields
        for table, table_obj in other_dataset.data.items():
            self._data[table] = table_obj.__class__(table, other_dataset.num_obs, dataset=self)
            self._data[table].copy_from(table_obj)
            self._fields.update({f: table for f in self._data[table].fields})

    def connect(self, table_from, table_to):
        """Connect two tables in the dataset

        TODO:
            This is currently not used.
        """
        self._data[table_from].connect_to(table_to)
        self._fields.update({f: table_from for f in self._data[table_from].fields})
        self._data[table_to].connect_from(table_from)
        self._fields.update({f: table_to for f in self._data[table_to].fields})

    def get_table(self, table):
        """Get a table of data from the dataset

        Note: This returns the actual data in the table, not a Table-object.

        Args:
            table:   String with name of table.

        Returns:
            Table data. The datatype will depend on the table.
        """
        return self._data[table].data

    def table_fields(self, table):
        """List names of fields in dataset belonging to a given table

        Args:
            table (String):  Name of table.

        Returns:
            List of strings with names of fields in table.
        """
        return sorted(f for f, t in self._fields.items() if t == table)

    def subset(self, idx):
        """Remove observations from all fields based on idx

        Modifies all fields in place, and updates the num_obs property.

        Args:
            idx:   Array of booleans with shape (num_obs, ). True means observation is kept, False means removed.
        """
        for table in self._data.values():
            table.subset(idx)
        self._num_obs = int(np.sum(idx))

    def extend(self, other_dset, meta_key=None):
        """Add observations from another dataset at the end of this dataset

        Note that this is quite strict in terms of which datasets can extend each other. They must have exactly the
        same tables. Tables containing several independent fields (e.g. float, text) may have different fields, in
        which case fields from both datasets are included. If a field is only defined in one dataset, it will get
        "empty" values for the observations in the dataset it is not defined.

        Args:
            other_dset (Dataset):  The other dataset.
            meta_key (str):        Dictionary key for introduction of an additional level in dictionary.
        """
        # Make sure tables are equal
        tbl_s = tuple([(k, type(v)) for k, v in sorted(self.data.items())])
        tbl_o = tuple([(k, type(v)) for k, v in sorted(other_dset.data.items())])
        if tbl_s != tbl_o:
            log.fatal(
                f"Dataset {self.description!r} can not be extended by {other_dset.description!r} "
                f"as their tables are different"
            )

        # Extend each table by calling extend on each of them
        for table_name, table in self._data.items():
            table.extend(other_dset.data[table_name])
        self._num_obs += other_dset.num_obs
        # TODO hjegei: Should it be done like that? If key already exists, should the values be merged together?
        if meta_key is None:
            self.meta.update(other_dset.meta)
        else:
            self.meta.setdefault(meta_key, dict()).update(other_dset.meta)

    def unique(self, field, **filters):
        """List all unique values of a given field

        TODO: Should this be implemented at the table level? Using np.unique? Seems like np.unique is faster

        Args:
            field:   String, fieldname.

        Returns:
            List of values. The value type depends on the field.
        """
        idx = [self.filter(**filters)] if filters else np.ones(self.num_obs, dtype=bool)
        if field in self._fields:
            return sorted(set(self[field][idx]))
        elif self.default_field_suffix and field + self.default_field_suffix in self._fields:
            return sorted(set(self[field + self.default_field_suffix][idx]))
        else:
            all_fields = [f for f in self._fields if f.startswith(field + "_")]
            return sorted(set().union(*[self[f][idx] for f in all_fields]))

    def filter(self, idx=None, **filters):
        """Filter observations

        Example:
            > dset.filter(station='NYALES20')
            array([True, True, False, ..., False, False, True], dtype=bool)

        Args:
            idx:      Optional, index with already filtered observations
            filters:  Key-value pairs representing field names and values.

        Returns:
            Array of booleans of shape (num_obs, ).
        """
        if idx is None:
            idx = np.ones(self.num_obs, dtype=bool)

        for field, filter_value in filters.items():
            if field in self._fields:
                table = self._fields[field]
                field_idx = self._data[table].filter_field(field, filter_value)
            elif self.default_field_suffix and field + self.default_field_suffix in self._fields:
                table = self._fields[field + self.default_field_suffix]
                field_idx = self._data[table].filter_field(field + self.default_field_suffix, filter_value)
            else:
                field_idx = np.zeros(self.num_obs, dtype=bool)
                for or_field in [f for f in self._fields if f.startswith(field + "_")]:
                    table = self._fields[or_field]
                    field_idx = np.logical_or(field_idx, self._data[table].filter_field(or_field, filter_value))
            idx = np.logical_and(idx, field_idx)
        return idx

    def first(self, field, **filters):
        """First value of field satisfying an optional filter

        """
        idx = np.where(self.filter(**filters))[0]
        if len(idx) == 0:
            return None

        return self[field][idx[0]]

    def for_each(self, key):
        previous_field_suffix = self.default_field_suffix
        suffixes = [f[len(key) :] for f in self.fields if f.startswith(key) and "." not in f]
        multipliers = [-1, 1] if len(suffixes) == 2 else [1] * len(suffixes)
        for multiplier, suffix in zip(multipliers, suffixes):
            self.default_field_suffix = suffix
            yield multiplier

        self.default_field_suffix = previous_field_suffix

    def values(self, *fields, filter=None):
        """Yield the values of fields for all observations

        Args:
            fields:   Strings representing field names.
            filter:   Dictionary with filter. See the filter-method.

        Returns:
            Iterator with values of the given fields.
        """
        if filter is None:
            idx = np.ones(self.num_obs, dtype=bool)
        else:
            idx = self.filter(**filter)

        for idx_in, *values in zip(idx, *[self[f] for f in fields]):
            if idx_in:
                yield values

    def unit(self, field):
        """Unit of values in a field in dataset

        Args:
            field (String):  The field name.

        Returns:
            String:  The unit used for data in the field.
        """
        table = self._data[self._fields[field]]
        return table.unit(field)

    def set_unit(self, field, unit):
        table = self._data[self._fields[field]]
        table.set_unit(field, unit)

    def plot_values(self, field):
        """Return values of a field in a form that can be plotted

        Args:
            field:   String, the field name.

        Returns:
            Numpy-array that can be plotted by for instance matplotlib.
        """
        table = self._data[self._fields[field]]
        return table.plot_values(field)

    def as_dict(self, use_plot_values=False, fields=None):
        """Return a representation of the dataset as a Dictionary

        Args:
            use_plot_values (Boolean):  Use the plot_values instead of regular values for each field.
            fields (List):              Field names that should be included, default is to include all fields.

        Returns:
            Dictionary: A representation of the dataset as a dictionary.
        """
        fields_dict = dict()
        for table in self._data.values():
            fields_dict.update(table.as_dict(use_plot_values=use_plot_values, fields=fields))

        return fields_dict

    def as_dataframe(self, index=None, use_plot_values=False, fields=None):
        """Return a representation of the dataset as a Pandas DataFrame

        Args:
            index (String / List):      Name of field to use as index. May also be a list of strings.
            use_plot_values (Boolean):  Use the plot_values instead of regular values for each field.
            fields (List):              Field names that should be included, default is to include all fields.

        Returns:
            DataFrame: A representation of the dataset as a Pandas DataFrame.
        """
        df = pd.DataFrame.from_dict(self.as_dict(use_plot_values=use_plot_values, fields=fields))
        if index is not None:
            df.set_index(index, drop=True, inplace=True)

        return df

    def rms(self, field, **filters):
        """Calculate root-mean-square of a field

        An optional filter can be given as keyword arguments, in which case the RMS is calculated only for observations
        satisying the filter.

        Args:
            field (String):  The field name.
            filters:         An optional filter, key-value pairs representing field names and values.

        Returns:
            Float:  The root-mean-square of a field.
        """
        return self.apply(lambda val: np.sqrt(np.mean(np.square(val))), field, **filters)

    def mean(self, field, **filters):
        """Calculate mean of a field

        An optional filter can be given as keyword arguments, in which case the mean is calculated only for
        observations satisying the filter.

        Args:
            field (String):  The field name.
            filters:         An optional filter, key-value pairs representing field names and values.

        Returns:
            Float:  The mean of a field.
        """
        return self.apply(np.mean, field, **filters)

    def apply(self, func, field, **filters):
        """Apply a function to a field

        An optional filter can be given as keyword arguments, in which case the function is applied only for
        observations satisying the filter.

        Example:
            > dset.apply(func=np.max, field='residual', satellite='G28')
            0.6523773

        Args:
            func (Function): The function to apply.
            field (String):  The field name.
            filters:         An optional filter, key-value pairs representing field names and values.

        Returns:
            Float:  The result of applying the function to the field.
        """
        values = self[field][self.filter(**filters)] if filters else self[field]
        return func(values)

    @property
    def rundate(self):
        """The rundate of the dataset

        Returns:
            Date:  The model run date.
        """
        return self._rundate

    @rundate.setter
    def rundate(self, value):
        """Set rundate of dataset and update corresponding variables

        Args:
            value (Date):  The model run date.
        """
        self._rundate = value
        self.vars["rundate"] = value.strftime(config.FMT_date)
        self.vars.update(**config.date_vars(value))

    @property
    def dataset_name(self):
        """The name of the dataset

        The full name of the dataset is `{dataset_name}/{dataset_id:04d}`.

        Returns:
            String: The name of the dataset
        """
        return self.name.split("/")[0]

    @property
    def dataset_id(self):
        """The id of the dataset

        The full name of the dataset is `{dataset_name}/{dataset_id:04d}`.

        Returns:
            Int: The id of the dataset
        """
        return int(self.name.split("/")[-1])

    @property
    def parameters(self):
        """The parameters used when initializing/loading the dataset

        Returns:
            Dict: The parameters used when initializing/loading the dataset
        """
        return dict(
            rundate=self.rundate,
            tech=self.vars["tech"],
            stage=self.vars["stage"],
            dataset_name=self.dataset_name,
            dataset_id=self.dataset_id,
            **self._kwargs,
        )

    @property
    def data(self):
        return self._data

    @property
    def default_field_suffix(self):
        return self._default_field_suffix

    @default_field_suffix.setter
    def default_field_suffix(self, suffix):
        if suffix:
            suffix = str(suffix)
            if not suffix.startswith("_"):
                suffix = "_" + suffix
        else:
            suffix = None
        self._default_field_suffix = suffix

    @property
    def fields(self):
        """List names of fields in dataset

        Returns:
            List of strings with names of fields.
        """
        return sorted(self._fields)

    @property
    def num_obs(self):
        """Number of observations (= rows) in dataset

        Returns:
            Int, number of observations in dataset.
        """
        return self._num_obs

    @num_obs.setter
    def num_obs(self, value):
        """Set number of observations in dataset

        This is part of the initialization of a dataset, and must be set before fields are added to the dataset. The
        number of observations can not be changed directly after fields are added to the dataset. In that case, one
        should use subset instead.

        Note that this means that it is impossible to add new observations to a dataset after fields are added to a
        dataset.

        Args:
            value:   Int, number of observations in dataset.
        """
        if self.fields:
            raise AttributeError(
                "Can not change number of observations after fields are added. " "Use 'subset' or 'extend' instead"
            )
        else:
            self._num_obs = value

    @property
    def tables(self):
        """List names of tables in dataset

        Returns:
            List of strings with names of tables.
        """
        return sorted(self._data)

    def add_to_meta(self, section, name, value):
        """Add meta information to the dataset.

        The meta-field is saved as a JSON-file so the values must be serializable to JSON (typically string, number,
        list or dictionary).

        Args:
            section (String):  Section inside meta-dict to store the value.
            name (String):     Name to associate value with.
            value:             Value to store in meta (datatype must be JSON-compatible).

        """
        meta_section = self.meta.setdefault(section, dict())
        meta_section[name] = value

    def add_event(self, timestamp, event_type, description):
        """Add information about an event to the dataset

        The events are stored in an `__events__`-field inside the meta-information of the dataset. Events can be
        retrieved using `get_events` or plotted in for instance There. The `event_type` can be used for grouping
        similar events together.

        Args:
            timestamp (Time):      Timestamp of event. Need not correspond to an observation epoch in the dataset.
            event_type (String):   Short description of the event type, for instance 'clock_break' or 'cycle_slip'.
            description (String):  More detailed description of the event.
        """
        events = self.meta.setdefault("__events__", dict())
        events.setdefault(event_type, list()).append((timestamp.utc.isot, description))

    def get_events(self, *event_types):
        """Retrieve events from the datasets

        One or several `event_types` can be specified, and only events of that type will be retrieved. If no event type
        is specified, all events are returned.

        TODO: Possible speed-up by creating one Time array instead of one object per event?
        TODO: Using set to remove duplicate events, how can these be handled so that duplicates are not written?

        Args:
            event_types (List):  Strings with event types like 'clock_break', 'cycle_slip', etc.

        Returns:
            List:  3-tuples describing each event, timestamp (Time), event_type (String), description (String).
        """
        all_events = self.meta.get("__events__", dict())
        if not event_types:
            event_types = all_events.keys()

        return sorted(
            set(
                [
                    (Time(e[0], format="isot", scale="utc"), k, e[1])
                    for k, v in all_events.items()
                    for e in v
                    if k in event_types
                ]
            )
        )

    def __getitem__(self, key):

        """Get a field from the dataset using dict-notation

        This method is called when a field is accessed using square brackets, e.g. data['station']. In general, we do a
        quick lookup in the _fields-dict and call __getitem__ on the correct table-object.

        There is also special handling of fields with a dot (.) in the field name.

        TODO:
            Explain the special handling of the dot!

        Args:
            key:   String, the name of the field.

        Returns:
            Field data. The datatype will depend on the table type.
        """
        # Add default field suffix if necessary
        if key not in self._fields and self.default_field_suffix is not None:
            key_list = key.split(".", maxsplit=1)
            key_list[0] += self.default_field_suffix
            key = ".".join(key_list)

        # Read field from proper table
        table = self._fields[key]
        if "." in key:
            key, attr = key.split(".", maxsplit=1)
            return getattr(self._data[table][key], attr)
        else:
            return self._data[table][key]

    def __delitem__(self, key):
        """Delete a field in the dataset

        This method is called by the del keyword when the field is specified using the square bracket notation, e.g.
        `del data['station']`. If the key is one of the fields we delete it from the field list together with any
        connected fields. If all fields in a table are deleted, we also delete the table from self._data.

        Note that we do not delete data from tables, which means that some "ghost data" might be saved by the tables
        (and written to disk), if not all fields in a table are deleted.

        TODO:
            Also delete actual data from tables.

        Args:
            key:   String, the name of the field.
        """
        if key in self._fields:
            for field in [f for f in self._fields if f == key or f.startswith(key + ".")]:
                table = self._fields[field]
                del self._fields[field]
                del self._data[table][field]
            for table in set(self._data) - set(self._fields.values()):
                del self._data[table]
        else:
            raise KeyError("{} '{}' has no field '{}'".format(type(self).__name__, self.dataset_name, key))

    def __getattr__(self, key):
        """Get a field from the dataset using the attribute-notation

        This method is called when a field is accessed using a dot, e.g. data.station. This is done by forwarding the
        call to __getitem__.

        Args:
            key:   String, the name of the field.

        Returns:
            Field data. The datatype will depend on the table type.
        """
        try:
            return self.__getitem__(key)
        except KeyError:
            error_msg = "'{}' object has no attribute '{}'".format(type(self).__name__, key)
            completions = [f for f in self.fields if f.startswith(key)]
            completion_msg = console.fill(
                ".\nMaybe one of these: " + ", ".join(completions),
                hanging=20,
                replace_whitespace=False,
                break_on_hyphens=False,
            )
            raise AttributeError(error_msg + (completion_msg if completions else "")) from None

    def __setattr__(self, key, value):
        """Make sure fields are not accidentally overwritten

        We override __setattr__ to intercept attempts at overwriting fields on the dataset, e.g. data.residual = obs -
        calc, because this will lead to subtle bugs where fields will seem like they are ok, but they will not update
        the underlying tables or be written to disk properly. Storing data to fields should typically be done with
        slice-notation instead: data.residual[:] = obs - calc.

        Technical note: Overriding __setattr__ comes with all kinds of consequences, since Python calls this method
        whenever an attribute is set. For this implementation to work it is necessary that _fields is the first
        attribute to be set on the Dataset-object (as is done in __init__). If _fields is not the first attribute to be
        set, the implicit getattr(self, '_fields') will fail with a KeyError and create an infinite loop. For the same
        reason, it is necessary to first call super().__setattr__ and later delete the key if it should not have been
        set.

        Args:
            key:   String, the name of the attribute.
            value: The new value of the attribute.
        """
        super().__setattr__(key, value)
        if key in self._fields:
            raise AttributeError(
                "Can't set attribute '{}' on '{}' object, use dset.{}[:] instead"
                "".format(key, type(self).__name__, key)
            )

    def __delattr__(self, key):
        """Delete a field in the dataset

        This method is called by the del keyword when a field is specified using the attribute notation, e.g.
        `del data.station`. If the key is one of the fields we delete it from the field list together with any
        connected fields. If all fields in a table are deleted, we also delete the table from self._data.

        Note that we do not delete data from tables, which means that some "ghost data" might be saved by the tables
        (and written to disk), if not all fields in a table are deleted.

        Args:
            key:   String, the name of the field.
        """
        try:
            self.__delitem__(key)
        except KeyError:
            raise AttributeError("'{}' object has no attribute '{}'".format(type(self).__name__, key)) from None

    def __dir__(self):
        """List all fields and attributes on the dataset

        Overridden so that fields are included in the list. This is especially useful when exploring the dataset
        interactively in an environment that supports tab-completion, e.g. ipython.

        Returns:
            List of strings, representing attributes and fields on the dataset.
        """
        return [f + ("(" if callable(getattr(self, f)) else "") for f in super().__dir__()] + self.fields

    def _ipython_key_completions_(self):
        """Tab completion for dict-style lookups in ipython

        See http://ipython.readthedocs.io/en/stable/config/integrating.html

        Returns:
            List: Strings with fieldnames in dataset.
        """
        return self.fields

    def __len__(self):
        """Length of dataset

        We define the length of the dataset as the number of observations in the dataset.

        Returns:
            Int: Number of observations in dataset.
        """
        return self.num_obs

    @property
    def description(self):
        """A one-line string describing the dataset.

        Returns:
            String: Description of the dataset.
        """
        return "{rundate} {tech}/{stage} {name}".format(
            rundate=self.rundate.strftime(config.FMT_date),
            tech=self.vars["tech"],
            stage=self.vars["stage"],
            name=self.name,
        )

    @property
    def repr(self):
        """A string representing the dataset

        The string should equal a python call that creates the Dataset. This is true if the dataset is written to disk.

        Note: This is implemented as the property `repr` instead of the special method `__repr__` as would be more
        conventional. This is a pragmatic choice to make it easier and more informational to work with datasets
        interactively.

        Returns:
            A string representing the dataset.
        """
        return "{cls}({rundate}, tech='{tech}', stage='{stage}', dataset_name='{name}', dataset_id={id})" "".format(
            cls=self.__class__.__name__,
            rundate=repr(self.rundate).split(".")[-1],
            tech=self.vars["tech"],
            stage=self.vars["stage"],
            name=self.dataset_name,
            id=int(self.dataset_id),
        )

    def __repr__(self):
        """A string describing the dataset

        The string describes the dataset and lists all fields. This string is mainly meant to be used when working
        interactively.

        Note: This is more typically implemented as `__str__` instead of `__repr__`. This is a pragmatic choice to make
        it easier and more informational to work with datasets. The usual `__repr__` string is available in the
        `repr`-property, while we use `__str__` to provide even more information.

        Returns:
            A string describing the dataset.
        """
        description = "{description}: {num_obs} obs\n".format(description=self.description, num_obs=self.num_obs)
        return description + console.fill("Fields: " + ", ".join(self.fields), hanging=8)

    def __str__(self):
        """A string describing the dataset

        The string describes the dataset, including listing all fields and their datatypes. This string is mainly
        meant to be read by humans.

        Returns:
            A string describing the dataset.
        """
        str_list = [
            "{}: {} obs, {} fields, {} tables".format(
                self.description, self.num_obs, len(self.fields), len(self.tables)
            )
        ]
        for _, table in sorted(self._data.items()):
            str_list.append(str(table))

        return "\n".join(str_list)

    @staticmethod
    def _add(datatype_cls):
        """Create a function that adds a field of a given datatype to the dataset

        Args:
            datatype_cls:  The datatype class. This should be a subclass of Table.

        Returns:
            Function that will add a field of the given datatype to the dataset.
        """

        def add_func(self, fieldnames, table=None, **kwargs):
            """This docstring is overwritten by the docstring in datatype_cls.add

            Args:
                fieldnames:  String or list of strings with names of fields to be added.
                table:       Name of table where fields are added.
            """
            if self.num_obs is None:
                raise InitializationError(
                    "Dataset must first be initialized by setting num_obs " "(number of observations)"
                )

            # Handle both one fieldname and a list of fieldnames by treating them as a list
            if isinstance(fieldnames, str):
                fieldnames = [fieldnames]

            for fieldname in fieldnames:
                if fieldname in self._fields:
                    raise FieldExistsError("The field '{}' already exists in Dataset".format(fieldname))
                if table is None:
                    table = datatype_cls.get_default_tablename(fieldname)
                if table not in self._data:
                    self._data[table] = datatype_cls(table, self.num_obs, dataset=self)

                # Add field to table and register field names
                self._data[table].add(fieldname, **kwargs)
                self._fields.update({f: table for f in self._data[table].fields})

        add_func.__doc__ = datatype_cls.add.__doc__
        return add_func

    @staticmethod
    def _read(datatype_cls):
        """Create a function that reads a table from file and adds it to the dataset

        Args:
            datatype_cls:  The datatype class. This should be a subclass of Table.

        Returns:
            Function that will read a table from file.
        """

        def read_func(self, table, json_data, hdf5_data):
            """Read data of type datatype_cls

            Args:
                table:     String, name of table where data are added.
                json_data: Dict, data read from JSON-file.
                hdf5_data: HDF5 dataset, data read from HDF5-file.
            """
            # Create table, read data and register fields
            self._data[table] = datatype_cls(table, self.num_obs, dataset=self)
            self._data[table].read(json_data, hdf5_data)
            self._fields.update({f: table for f in self._data[table].fields})

        return read_func
