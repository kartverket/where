"""A Dataset direction field"""

# Third party imports
import numpy as np

# Midgard imports
from midgard.data.fieldtypes._fieldtype import FieldType
from midgard.dev import exceptions
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where.data.direction import Direction
from where.data._direction import DirectionArray


@plugins.register
class DirectionField(FieldType):

    _subfields = DirectionArray.fieldnames()
    _plotfields = DirectionArray.plot_fields()
    _factory = staticmethod(Direction)

    def _post_init(self, val, **dir_args):
        """Initialize float field"""
        if isinstance(val, DirectionArray):
            data = val
        else:
            data = self._factory(val, **dir_args)

        # Check that unit is not given, overwrite with direction units
        if self._unit is not None and self._unit != data.unit():
            raise exceptions.InitializationError("Parameter 'unit' should not be specified for directions")
        self._unit = data.unit()

        # Check that the correct number of observations are given
        if len(data) != self.num_obs:
            raise ValueError(f"{self.name!r} initialized with {len(data)} values, expected {self.num_obs}")

        # Check that the correct number of columns are given
        if data.ndim != 2:
            raise ValueError(f"{self.name!r} initialized with {data.ndim} columns, expected 2 (ra, dec)")

        # Store the data as a TimeArray
        self.data = data

    def plot_values(self, field=None) -> np.array:
        """Return values of the field in a form that can be plotted"""
        if not field:
            return self.data.val

        values = getattr(self.data, field)
        if isinstance(values, DirectionArray):
            return values.val
        else:
            return values

    def _prepend_empty(self, num_obs, memo):
        empty_shape = (num_obs, *self.data.shape[1:])
        empty = DirectionArray(np.full(empty_shape, np.nan))

        empty_id = id(empty)
        self.data = DirectionArray.insert(self.data, 0, empty, memo)
        memo.pop(empty_id, None)

    def _append_empty(self, num_obs, memo):
        empty_shape = (num_obs, *self.data.shape[1:])
        empty = DirectionArray(np.full(empty_shape, np.nan))

        empty_id = id(empty)
        self.data = DirectionArray.insert(self.data, self.num_obs, empty, memo)
        memo.pop(empty_id, None)

    def _subset(self, idx, memo):
        self.data = self.data.subset(idx, memo)

    def _extend(self, other_field, memo) -> None:
        """Add observations from another field"""
        if other_field.data.ndim != self.data.ndim:
            raise ValueError(
                f"Field '{self.name}' cannot be extended. Dimensions must be equal. ({other_field.data.ndim} != {self.data.ndim})"
            )

        try:
            factors = [Unit(from_unit, to_unit) for from_unit, to_unit in zip(other_field._unit, self._unit)]
        except exceptions.UnitError:
            raise exceptions.UnitError(
                f"Cannot extend field '{self.name}'. {other_field._unit} cannot be converted to {self._unit}"
            )
        except TypeError:
            if self._unit == other_field._unit == None:
                factors = 1
            else:
                raise exceptions.UnitError(
                    f"Cannot extend field '{self.name}'. {other_field._unit} cannot be converted to {self._unit}"
                )

        self.data = DirectionArray.insert(self.data, self.num_obs, other_field.data * factors, memo)

    @classmethod
    def _read(cls, h5_group, memo) -> "DirectionField":
        """Read a DirectionField from a HDF5 data source"""
        name = h5_group.attrs["fieldname"]
        if name in memo:
            direction = memo[name]
        else:
            direction = DirectionArray._read(h5_group, memo)

        return cls(num_obs=len(direction), name=name, val=direction)

    def _write(self, h5_group, memo) -> None:
        """Write a DirectionField to a HDF5 data source"""
        self.data._write(h5_group, memo)
