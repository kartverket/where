"""Framework for reading and working with Celestial Reference Frames

Description:
------------

Each type of reference frame should be defined as a class in a separate .py-file. The class inside the .py-file that
represents the reference frame need to be decorated with the :func:`~midgard.dev.plugins.register` decorator as follows::

    from midgard.dev import plugins

    @plugins.register
    class Icrf2(crf.CrfFactory):
        ...

The reference frame class will be responsible for reading the necessary data using the proper parsers and for providing
radio source information. It must inherit from apriori.crf.CrfFactory.

To obtain a CRF, call :func:`apriori.get('crf', ...)`. Which reference frames to use can be specified explicitly. If it
is not, the reference frames are looked up in the config files. See :func:`get_crf` for more information.




"""
# Where imports
from where.lib import config
from midgard.dev import plugins

# Make Crf-objects available for convenience
from where.apriori.crf._crf import Crf, CrfFactory, CrfSource  # noqa


def get_crf(time, celestial_reference_frames=None):
    """Get a celestial reference frame

    The specification of the reference frames takes the form `crf1, crf2, ...` where `crf1` etc is matched with a
    filename in this crf-directory. The version directive can be interpreted differently by each reference frame. Here
    is one example of some typical reference frame specifications::

        celestial_reference_frames = icrf2                 # The ICRF version 2
        celestial_reference_frames = vlbi_obs              # The radio source info from the observation file
        celestial_reference_frames = vascc                 # The CRF used in VASCC

    The reference frames are used in the order given. For instance if::

        celestial_reference_frames = Icrf2, vlbi_obs

    then all radio sources will be looked up in the Icrf2-frame. If a radio source is not found there, it will be
    looked up in the vlbi_obs-reference frame instead.

    Args:
        time (Time):                Time epochs for which to calculate the reference frame.
        celestial_reference_frames (String):  Optional specification of which reference frames to use (see above).

    Returns:
        Crf: Reference frame object.

    """
    celestial_reference_frames = config.tech.get("celestial_reference_frames", celestial_reference_frames).str
    return Crf(time=time, celestial_reference_frames=celestial_reference_frames)


def get_crf_factory(time, celestial_reference_frame):
    """Get a factory for a given celestial reference frame

    The factory knows how to create RadioSource objects for a given reference frame, for instance `icrf2`.

    Args:
        celestial_reference_frame (String):  Specification of which reference frame to use (see `get_crf`).

    Returns:
        CrfFactory:  Factory that knows how to create RadioSource objects.
    """
    name, _, catalog = celestial_reference_frame.partition(":")
    kwargs = dict(catalog=catalog) if catalog else dict()
    return plugins.call(package_name=__name__, plugin_name=name, time=time, **kwargs)
