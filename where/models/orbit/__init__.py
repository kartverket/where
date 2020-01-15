"""Framework for calculating satellite orbit models

Description:
------------

Each orbit model should be defined in a separate .py-file. The function inside the .py-file that should be called need
to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    def gravity_earth(...):
        ...

The decorated function will be called with several parameters::

    TODO: Document parameters
    TODO: Plugins does not work with cython

References:
-----------

[1] Montenbruck, Oliver and Gill, Eberhard: Satellite Orbits, Springer Verlag, 2000.




"""

from where.models.orbit._orbit import calculate as calculate_orbit, update_orbit  # noqa
