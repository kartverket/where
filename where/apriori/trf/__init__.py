"""Framework for reading and working with Terrestrial Reference Frames

Description:
------------

Each type of reference frame should be defined as a class in a separate .py-file. The class inside the .py-file that
represents the reference frame need to be decorated with the :func:`~where.lib.plugins.register` decorator as follows::

    from where.lib import plugins

    @plugins.register
    class Itrf(trf.TrfFactory):
        ...

The reference frame class will be responsible reading the necessary data using the proper parsers and for providing the
velocity model. It must inherit from apriori.trf.TrfFactory.

To obtain a TRF, call :func:`apriori.get('trf', time=time, ...)` with the requested time epochs. Which
reference frames to use can be specified explicitly. If it is not, the reference frames are looked up in the config
files. See :func:`get_trf` for more information.




"""
# Where imports
from where.lib import config
from where.lib import plugins

# Make Trf-objects available for convenience
from where.apriori.trf._trf import Trf, TrfFactory, TrfSite  # noqa


def get_trf(time, reference_frames=None):
    """Get a reference frame for the given time epochs

    The specification of the reference frames takes the form `trf1:version, trf2:version` where `trf` is matched with a
    filename in this trf-directory. The version directive can be interpreted differently by each reference frame. Here
    is one example of some typical reference frame specifications::

        reference_frames = itrf:2014             # The ITRF 2014 reference frame (sinex file)
        reference_frames = itrf:2014_ssc         # The ITRF 2014 reference frame (ssc file)
        reference_frames = itrf:2008_ssc         # The ITRF 2008 reference frame (ssc file)

    The reference frames are used in the order given. For instance if::

        reference_frames = itrf:2014, vlbi_obs

    then all sites will be looked up in the itrf:2014-frame. If a site is not found there, it will be looked up in the
    vlbi_obs-reference frame instead.

    Args:
        time (Time):                Time epochs for which to calculate the reference frame.
        reference_frames (String):  Optional specification of which reference frames to use (see above).

    Returns:
        Trf: Reference frame object.
    """
    reference_frames = config.tech.get("reference_frames", reference_frames).str
    return Trf(time=time, reference_frames=reference_frames)


def get_trf_factory(time, reference_frame):
    """Get a factory for a given reference frame

    The factory knows how to create TrfSite objects for a given reference frame, for instance `itrf:2014`.

    Args:
        time (Time):               Time epochs for which to calculate the reference frame.
        reference_frame (String):  Specification of which reference frame to use (see `get_trf`).

    Returns:
        TrfFactory:  Factory that knows how to create TrfSite objects.
    """
    name, _, version = reference_frame.partition(":")
    kwargs = dict(version=version) if version else dict()
    return plugins.call_one(package_name=__name__, plugin_name=name, time=time, **kwargs)
