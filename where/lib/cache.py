"""Where library module for caching

Example:
    from where.lib import cache


Description:

asdf.

Ref Python Cookbook, 3rd ed. Recipe 8.10.


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $
"""
# Standard library imports
import functools


_DEPENDENT_PROPERTIES = dict()


class property():

    def __init__(self, fget):
        self.fget = fget
        functools.update_wrapper(self, fget)

    def __get__(self, instance, cls):
        if instance is None:
            return self
        else:
            value = self.fget(instance)
            setattr(instance, self.fget.__name__, value)
            return value


class register_dependencies(type):

    _dependencies = list()

    def __call__(cls, fget):
        *container, name = fget.__qualname__.split(".")
        for dependency in cls._dependencies:
            _DEPENDENT_PROPERTIES.setdefault(dependency, set()).add((".".join(container), name))
        return super().__call__(fget)

    def __getattr__(cls, key):
        cls._dependencies.append(key)
        return cls


class dependent_property(property, metaclass=register_dependencies):
    pass


def forget_dependent_values(obj, *dependencies):
    for dependency in dependencies:
        for container, prop in _DEPENDENT_PROPERTIES.get(dependency, ()):
            if container == obj.__class__.__qualname__:
                try:
                    delattr(obj, prop)
                except AttributeError:
                    pass


def function(func):
    """Cache a given function call

    Uses the lru_cache (Least Recently Used) implementation from the functools standard library. Wrapped here to have a
    consistent lib.cache-module and possibly set the parameters of the functools.lru_cache.
    """
    return functools.lru_cache()(func)
