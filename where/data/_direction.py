"""Array with direction values

See https://docs.scipy.org/doc/numpy/user/basics.subclassing.html for information about subclassing Numpy arrays.

DirectionArray is a regular Numpy array with two columns, right ascension and direction
"""

# Standard library imports
from typing import Any, Callable, Dict, List, Tuple
import sys
import weakref
from functools import lru_cache


# Third party imports
import numpy as np

# Midgard imports
from midgard.dev import exceptions

_ATTRIBUTES: List[str] = list()  # Populated by register_attribute()
_FIELDS: List[str] = list()  # Populated by register_field()
_UNITS: Dict[str, str] = dict()  # Populated by register_field()
_SYSTEMS: Dict[str, Callable] = dict()  # Populated by register_system()
_CONVERSIONS: Dict[Tuple[str, str], Callable] = dict()  # Populated by register_system()
_CONVERSION_HOPS: Dict[Tuple[str, str], List[str]] = dict()  # Cache for to_system()


def register_attribute(name: str) -> None:
    """Function used to register new attributes on position arrays

    The registered attributes will be available as attributes on PositionArray
    and its subclasses. In addition, each attribute can be given as a parameter
    when creating a PositionArray.

    The reason for using this register-function instead of a regular attribute
    is to allow additional attributes to be added on all position systems.

    Args:
        name:     Name of attribute
        attr_cls: Class of attribute

    """
    _ATTRIBUTES.append(name)


def register_system(
    convert_to: Dict[str, Callable] = None, convert_from: Dict[str, Callable] = None
) -> Callable[[Callable], Callable]:
    """Decorator used to register new position systems

    The system name is read from the .system attribute of the Position class.

    Args:
        convert_to:    Functions used to convert to other systems.
        convert_from:  Functions used to convert from other systems.

    Returns:
        Decorator registering system.
    """

    def wrapper(cls: Callable) -> Callable:
        name = cls.system
        _SYSTEMS[name] = cls

        if convert_to:
            for to_system, converter in convert_to.items():
                _CONVERSIONS[(name, to_system)] = converter
        if convert_from:
            for from_system, converter in convert_from.items():
                _CONVERSIONS[(from_system, name)] = converter
        return cls

    return wrapper


def register_field(units: List[str]) -> Callable:
    """Decorator used to register fields and their units

    """

    def wrapper(func: Callable) -> Callable:
        field = func.__name__
        _FIELDS.append(field)
        _UNITS[field] = units
        return func

    return wrapper


def _find_conversion_hops(hop: Tuple[str, str]) -> List[Tuple[str, str]]:
    """Calculate the hops needed to convert between systems using breadth first search"""
    start_sys, target_sys = hop
    queue = [(start_sys, [])]
    visited = set()

    while queue:
        from_sys, hops = queue.pop(0)
        for to_sys in [t for f, t in _CONVERSIONS if f == from_sys]:
            one_hop = (from_sys, to_sys)
            if to_sys == target_sys:
                return hops + [one_hop]
            if one_hop not in visited:
                visited.add(one_hop)
                queue.append((to_sys, hops + [one_hop]))

    raise exceptions.UnknownConversionError(f"Can't convert PositionArray from {start_sys!r} to {target_sys!r}")


class DirectionArray(np.ndarray):

    type = "direction"
    _units = ("unitless", "unitless", "unitless")

    def __new__(cls, val, **dir_args):
        """Create a new DirectionArray"""
        obj = np.asarray(val, dtype=float, order="C").view(cls)
        obj.system = cls.system

        for attr in cls._attributes():
            setattr(obj, attr, dir_args.get(attr))
        return obj

    def __array_finalize__(self, obj):
        """Called automatically when a new DirectionArray is created"""
        self._cache = dict()
        self._dependent_objs = list()

        if obj is None:
            return

        # Validate shape
        num_columns = len(self.column_names)
        if self.shape[-1] != num_columns:
            column_names = ", ".join(self.column_names)
            raise ValueError(f"{type(self).__name__!r} requires {num_columns} columns: {column_names}")

        if self.ndim > 2:
            raise ValueError(f"{type(self).__name__!r} must be a 1- or 2-dimensional array with {num_columns} columns")

        # Copy attributes from the original object
        self.system = getattr(obj, "system", None)
        for attr in self._attributes():
            attr_sliced = getattr(obj, f"_{attr}_sliced", None)
            if attr_sliced is not None:
                setattr(self, attr, attr_sliced)
            else:
                setattr(self, attr, getattr(obj, attr, None))

    @staticmethod
    def create(val: np.ndarray, system: str, **dir_args: Any) -> "DirectionArray":
        """Factory for creating DirectionArrays for different systems

        See each direction class for exact optional parameters.

        Args:
            val:       Array of position values.
            system:    Name of direction system.
            dir_args:  Additional arguments used to create the DirectionArray.

        Returns:
            Array with positions in the given system.
        """
        if system not in _SYSTEMS:
            systems = ", ".join(_SYSTEMS)
            raise exceptions.UnknownSystemError(f"System {system!r} unknown. Use one of {systems}")

        return _SYSTEMS[system](val, **dir_args)

    def fieldnames(self):
        """Return list of valid attributes for this object"""
        # Pick one element to avoid doing calculations on a large array 
        obj = self if len(self) == 1 else self[0]

        systems_and_columns = []
        for system in obj._systems():
            try:
                _find_conversion_hops((obj.system, system))
                # Add systems
                systems_and_columns.append(system)
                for column in _SYSTEMS[system].column_names:
                    # Add system columns
                    systems_and_columns.append(f"{system}.{column}")
            except exceptions.UnknownConversionError:
                pass  # Skip systems that cannot be converted to

        useable_fields = []
        for field in self._fields():
            try:
                getattr(obj, field)
                useable_fields.append(field)
            except exceptions.InitializationError:
                # This field cannot be computed. Skipping
                pass

        return self._attributes() + systems_and_columns + useable_fields

    def plot_fields(self):
        """Return list of plottable attributes for this object"""
        valid_fields = set(self.fieldnames())
        return list(valid_fields - set(self._attributes()))

    @classmethod
    def _attributes(cls):
        return _ATTRIBUTES

    @classmethod
    def _fields(cls):
        return _FIELDS

    @classmethod
    def _systems(cls):
        return list(_SYSTEMS.keys())

    @classmethod
    def _system_columns(cls):
        return [f"{s}.{c}" for s, sc in _SYSTEMS.items() for c in sc.column_names]

    def to_system(self, system: str) -> "PosDeltaBase":
        """Convert to a different system

        Args:
            system:  Name of new system.

        Returns:
            PosDeltaBase representing the same positions or position deltas in the new system.
        """
        # Don't convert if not necessary
        if system == self.system:
            return self

        # Raise error for unknown systems
        if system not in _SYSTEMS:
            systems = ", ".join(_SYSTEMS)
            raise exceptions.UnknownSystemError(f"System {system!r} unknown. Use one of {systems}")

        if system in self._cache:
            return self._cache[system]

        # Convert to new system
        hop = (self.system, system)
        if hop in _CONVERSIONS:
            self._cache[system] = _SYSTEMS[system].convert_to(self, _CONVERSIONS[hop])
            return self._cache[system]

        if hop not in _CONVERSION_HOPS:
            _CONVERSION_HOPS[self.cls_name][hop] = _find_conversion_hops(hop)

        val = self
        for one_hop in _CONVERSION_HOPS[hop]:
            val = _SYSTEMS[one_hop[-1]].convert_to(val, _CONVERSIONS[one_hop])
            self._cache[one_hop[-1]] = val
        return val

    @classmethod
    def convert_to(cls, direction: "DirectionArray", converter: Callable) -> "DirectionArray":
        """Convert the direction to the type of this class

        Applies the converter function that is provides and copies all registered attributes to the new position 
        """
        attrs = {a: getattr(direction, a, None) for a in cls._attributes()}
        return _SYSTEMS[cls.system](converter(direction), **attrs)

    @property
    def is_transposed(self):
        # Because we forced order == "C" on creation
        return self.flags.f_contiguous

    @property
    def pos(self):
        """PositionArray.other needs to have a pos attribute"""
        return self

    @property
    def val(self):
        return np.asarray(self)

    @classmethod
    def unit(cls, field: str = "") -> Tuple[str, ...]:
        """Unit of field"""
        mainfield, _, subfield = field.partition(".")

        # Units of direction array
        if not field:
            return cls._units
        # Unit of columns in direction array
        elif field in cls.column_names:
            return (cls._units[cls.column_names.index(field)],)
        # Units of properties
        elif field in _UNITS:
            return _UNITS[field]
        # Units of systems
        elif mainfield in _SYSTEMS:
            return _SYSTEMS[mainfield].unit(subfield)
        else:
            raise exceptions.FieldDoesNotExistError(f"Field {mainfield!r} does not exist") from None

    @property
    @lru_cache()
    @register_field(units=("radians",))
    def right_ascension(self):
        """Right ascension"""
        if "right_ascension" not in self._cache:
            ra = np.arctan2(self.y, self.x)
            if self.ndim == 1:
                if ra < 0:
                    ra = ra + 2 * np.pi
            else:
                ra[ra < 0] += 2 * np.pi
            self._cache["right_ascension"] = ra
        return self._cache["right_ascension"]

    @property
    @lru_cache()
    @register_field(units=("radians",))
    def declination(self):
        """Declination"""
        if "declination" not in self._cache:
            self._cache["declination"] = np.pi / 2 - np.arccos(self.z)
        return self._cache["declination"]

    @property
    @lru_cache()
    @register_field(units=("unitless", "unitless", "unitless"))
    def dsrc_dra(self):
        """Derivative of radio source unit vector with regards to right ascension"""
        if "dsrc_dra" not in self._cache:
            cos_dec = np.cos(self.declination)
            cos_ra = np.cos(self.right_ascension)
            sin_ra = np.sin(self.right_ascension)
            zero = np.zeros_like(cos_dec)
            self._cache["dsrc_dra"] = np.array([-cos_dec * sin_ra, cos_dec * cos_ra, zero]).T
        return self._cache["dsrc_dra"]

    @property
    @lru_cache()
    @register_field(units=("unitless", "unitless", "unitless"))
    def dsrc_ddec(self):
        """Derivative of radio source unit vector with regards to declination"""
        if "dsrc_ddec" not in self._cache:
            cos_dec = np.cos(self.declination)
            cos_ra = np.cos(self.right_ascension)
            sin_dec = np.sin(self.declination)
            sin_ra = np.sin(self.right_ascension)
            self._cache["dsrc_ddec"] = np.array([-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec]).T
        return self._cache["dsrc_ddec"]

    @property
    @register_field(units=("unitless", "unitless", "unitless"))
    def unit_vector(self):
        """Unit vector"""
        return np.asarray(self)

    def direction_from(self, _):
        """Calcualte direction vector from position ignoring the motion of the Earth"""
        # Ignore station position
        return self.unit_vector

    def __hash__(self):
        return hash(self.tobytes())

    def __eq__(self, other):
        return self.data.tobytes() == other.data.tobytes()

    def __deepcopy__(self, memo):
        new_sigma_array = self.__class__(np.asarray(self))
        memo[id(self)] = new_sigma_array
        return new_sigma_array

    def __getattr__(self, key):
        """Get attributes with dot notation

        Add systems and column names to attributes on Position and PostionDelta arrays.

        Args:
            key:  Name of attribute.

        Returns:
            Value of attribute.
        """
        mainfield, _, subfield = key.partition(".")

        if subfield:
            return getattr(getattr(self, mainfield), subfield)

        # Convert to a different system
        if key in _SYSTEMS:
            return self.to_system(key)

        # Return one column as a regular numpy array
        elif key in self.column_names:
            idx = self.column_names.index(key)
            if self.ndim == 1:
                return self.val[idx]
            else:
                if self.is_transposed:
                    return self.val[idx, :]
                else:
                    return self.val[:, idx]

        # Raise error for unknown attributes
        else:
            raise AttributeError(f"{type(self).__name__!r} has no attribute {key!r}") from None

    def __getitem__(self, item):
        """Update attributes with correct shape, used by __array_finalize__"""

        # Get column
        if isinstance(item, int) and self.is_transposed:
            return getattr(self, self.column_names[item])

        from_super = super().__getitem__(item)

        # Get row
        if isinstance(item, int):
            dir_args = {}
            for attr in self._attributes():
                orig_value = getattr(self, attr, None)
                if orig_value is not None:
                    sliced_value = orig_value[item]
                    sliced_attr = f"_{attr}_sliced"
                    setattr(self, sliced_attr, sliced_value)
                    dir_args[attr] = sliced_value
            return self.__class__(from_super, **dir_args)

        return from_super

    @classmethod
    def insert(cls, a, pos, b, memo):
        """ Insert b into a at position pos"""
        id_a = id(a)
        if id_a in memo:
            return memo[id_a][-1]

        id_b = id(b)
        if id_b in memo:
            return memo[id_b][-1]

        val = np.insert(np.asarray(a), pos, np.asarray(b))
        new_direction = cls(val)
        memo[id_a] = (a, new_direction)
        memo[id_b] = (b, new_direction)
        return new_direction

    def subset(self, idx, memo):
        """Create a subset """
        old_id = id(self)
        if old_id in memo:
            return memo[old_id]

        val = np.asarray(self)[idx]
        dir_args = dict()

        for attr_name in self._attributes():
            attr = getattr(self, attr_name, None)
            if attr is None:
                continue

            old_id_attr = id(attr)
            if old_id_attr in memo:
                dir_args[attr_name] = memo[old_id_attr]
            else:
                dir_args[attr_name] = attr.subset(idx, memo)
                memo[old_id_attr] = dir_args[attr_name]

        new_pos = _SYSTEMS[self.system](val, **dir_args)
        memo[old_id] = new_pos
        return new_pos

    def __matmul__(self, _):
        """self @ _"""
        return NotImplemented

    def __rmatmul__(self, _):
        """_ @ self"""
        return NotImplemented

    def __imatmul__(self, _):
        """self @= _"""
        return NotImplemented

    def __mul__(self, _):
        """self * other """
        return NotImplemented

    def __rmul__(self, other):
        """other * self"""
        return NotImplemented

    def __imul__(self, _):
        """self *= _"""
        return NotImplemented

    def __truediv__(self, _):
        """self / _"""
        return NotImplemented

    def __rtruediv__(self, _):
        """_ / self"""
        return NotImplemented

    def __itruediv__(self, _):
        """self /= _"""
        return NotImplemented

    def __floordiv__(self, _):
        """self // _"""
        return NotImplemented

    def __rfloordiv__(self, _):
        """ _ // self"""
        return NotImplemented

    def __ifloordiv__(self, _):
        """self //= _"""
        return NotImplemented

    def __pow__(self, _):
        """ self ** _"""
        return NotImplemented

    def __rpow(self, _):
        """ _ ** self """
        return NotImplemented

    def __ipow__(self, _):
        """ self **= _"""
        return NotImplemented

    def __setitem__(self, key, item):
        self.clear_cache()  # Clear cache when any elements change
        if self._dependent_objs:
            for o in self._dependent_objs:
                o().clear_cache()  # Clear cache of dependent obj
        return super().__setitem__(key, item)

    def __setattr__(self, key, value):
        self.clear_cache()
        if key in self._attributes():
            prev_attr_value = getattr(self, key, None)
            if prev_attr_value is not None:
                prev_attr_value.remove_dependency(self)
            if value is not None:
                try:
                    value.add_dependency(self)
                except AttributeError:
                    pass
        return super().__setattr__(key, value)

    def clear_cache(self):
        for k, v in getattr(self, "_cache", {}).items():
            if k in self._systems():
                try:
                    v.other.remove_dependency(v)
                except AttributeError:
                    pass
        super().__setattr__("_cache", dict())

    def add_dependency(self, dependency):
        self._dependent_objs.append(weakref.ref(dependency))
        # Clean out references that have been garbage collected

        dead_refs = [obj for obj in self._dependent_objs if obj() is None]
        for element in dead_refs:
            self._dependent_objs.remove(element)

    def remove_dependency(self, dependency):
        try:
            idx = next(idx for idx, obj in enumerate(self._dependent_objs) if id(obj()) == id(dependency))
            del self._dependent_objs[idx]
        except StopIteration:
            pass

    @classmethod
    def _read(cls, h5_group, memo):
        system = h5_group.attrs["system"]

        dir_args = {}
        for a in cls._attributes():
            if a in h5_group.attrs:
                # attribute is a reference to the data of another field
                fieldname = h5_group.attrs[a]
                if fieldname in memo:
                    # the other field has already been read
                    dir_args.update({a: memo[fieldname]})
                else:
                    # the other field has not been read yet
                    attr_group = h5_group.parent[fieldname]
                    cls_module, _, cls_name = attr_group.attrs["__class__"].rpartition(".")
                    attr_cls = getattr(sys.modules[cls_module], cls_name)
                    arg = attr_cls._read(attr_group, memo)
                    dir_args.update({a: arg})
                    memo[fieldname] = arg
            elif a in h5_group and isinstance(h5_group[a], type(h5_group)):
                # attribute is a part of this group and is not in a separate field
                cls_module, _, cls_name = h5_group[a].attrs["__class__"].rpartition(".")
                attr_cls = getattr(sys.modules[cls_module], cls_name)
                arg = attr_cls._read(h5_group[a], memo)
                dir_args.update({a: arg})
                memo[f"{h5_group.attrs['fieldname']}.{a}"] = arg

        val = h5_group[h5_group.attrs["fieldname"]][...]

        direction = cls.create(val, system=system, **dir_args)
        memo[f"{h5_group.attrs['fieldname']}"] = direction
        return direction

    def _write(self, h5_group, memo):
        h5_group.attrs["system"] = self.system
        h5_field = h5_group.create_dataset(h5_group.attrs["fieldname"], self.shape, dtype=self.dtype)
        h5_field[...] = self.val

        for a in self._attributes():
            attr = getattr(self, a, None)
            if attr is None:
                continue

            if id(attr) in memo:
                # attribute is a reference to the data of another field
                h5_group.attrs[a] = memo[id(attr)]
            else:
                # attribute is stored as part of this in PositionArray
                h5_sub_group = h5_group.create_group(a)
                h5_sub_group.attrs["fieldname"] = a
                h5_sub_group.attrs["__class__"] = f"{attr.__class__.__module__}.{attr.__class__.__name__}"
                memo[id(attr)] = f"{h5_group.attrs['fieldname']}.{a}"
                attr._write(h5_sub_group, memo)  # Potential recursive call
