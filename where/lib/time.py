"""Where library module for working with times and dates

Example:
    from where.lib import time
    t = time.Time([datetime(2016, 1, 18, h, 0) for h in range(24)], scale='utc')

Description:

This module contains functions and classes for working with times and dates in Where. Most important is the Time-class
which is a subclass of the Astropy Time-class, and behaves mostly identically.

References:
[1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010).
    http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

[2] Astropy Source Code
    https://github.com/astropy/astropy


$Revision: 15123 $
$Date: 2018-05-16 18:00:58 +0200 (Wed, 16 May 2018) $
$LastChangedBy: dahmic $

"""
# Standard library imports
from datetime import datetime

# External library imports
import astropy.time
import numpy as np

# Where imports
from where.lib import cache
from where.lib import constant
from where.ext import iers_2010 as iers
from where.ext import sofa_wrapper as sofa

# Make astropy.time.TimeDelta available
from astropy.time import TimeDelta  # noqa


class Time(astropy.time.Time):
    """Represent and manipulate times and dates in Where

    The implementation is based off the Astropy Time-class. The differences are that we set the time scale deltas ut1 -
    utc and tdb - tt automatically, we add 'gps' as a time scale and there are some more time formats available,
    including jd_int, jd_frac, mjd_int, mjd_frac, gpsweek, gpsday and gpssec (technically, these are not proper
    formats, but properties. This means that they can not be used for input). The old 'gps'-format has been renamed to
    'gps_1980'.

    Available time scales:

    ======  ===================================
    Scale 	Description
    ======  ===================================
    gps     GPS Time (GPS)
    tai     International Atomic Time (TAI)
    tcb     Barycentric Coordinate Time (TCB)
    tcg     Geocentric Coordinate Time (TCG)
    tdb     Barycentric Dynamical Time (TDB)
    tt      Terrestrial Time (TT)
    ut1     Universal Time (UT1)
    utc     Coordinated Universal Time (UTC)
    ======  ===================================

    Available Astropy time formats:

    ============    ===================================
    Format          Example argument
    ============    ===================================
    byear           1950.0
    byear_str       'B1950.0'
    cxcsec          63072064.184
    datetime        datetime(2000, 1, 2, 12, 0, 0)
    decimalyear     2000.45
    fits            '2000-01-01T00:00:00.000(TAI)'
    gps_1980        630720013.0
    iso             '2000-01-01 00:00:00.000'
    isot            '2000-01-01T00:00:00.000'
    jd              2451544.5
    jyear           2000.0
    jyear_str       'J2000.0'
    mjd             51544.0
    plot_date       730120.0003703703
    unix            946684800.0
    yday            2000:001:00:00:00.000
    ============    ===================================

    Additional Where time formats:

    ============  =========  ========================================
    Format        Example    Description
    ============  =========  ========================================
    jd_frac       0.5475     Fractional part of Julian Day
    jd_int        2451544.0  Integer part of Julian Day
    mjd_frac      0.0475     Frational part of Modified Julian Day
    mjd_int       51544.0    Integer part of Modified Julian Day
    gpsday        6          Day of GPS week
    gpssec        432000.0   Seconds of GPS week
    gpsweek       1042.0     GPS week
    ============  =========  ========================================

    """

    @cache.property
    def text_repr(self):
        """A text representation of all time epochs

        This is used when comparing two Time objects using either `__hash__` or `__eq__`.

        Returns:
            String: Text representation of all time epochs.
        """
        return self.scale + (self.isot if self.isscalar else "".join(self.isot.tolist()))

    @cache.property
    def itrs2gcrs(self):
        """The ITRS to GCRS transformation matrix for the given times

        See ext/sofa_wrapper.py for details about the implementation.

        Returns:
            Numpy-float array with transformation matrices, shape is (len(self), 3, 3).
        """
        return sofa.Q(self) @ sofa.R(self) @ sofa.W(self)

    @cache.property
    def gcrs2itrs(self):
        """The GCRS to ITRS transformation matrix for the given times

        See ext/sofa_wrapper.py for details about the implementation.

        Returns:
            Numpy-float array with transformation matrices, shape is (len(self), 3, 3).
        """
        if self.isscalar:
            return np.transpose(self.itrs2gcrs)
        else:
            return np.transpose(self.itrs2gcrs, axes=[0, 2, 1])

    @cache.property
    def itrs2gcrs_dot(self):
        """The derivative of ITRS to GCRS transformation matrix for the given times

        See ext/sofa_wrapper.py for details about the implementation.

        Returns:
            Numpy-float array with transformation matrices, shape is (len(self), 3, 3).
        """
        return sofa.Q(self) @ sofa.dR_dut1(self) @ sofa.W(self)

    @cache.property
    def gcrs2itrs_dot(self):
        """The derivative of GCRS to ITRS transformation matrix for the given times

        See ext/sofa_wrapper.py for details about the implementation.

        Returns:
            Numpy-float array with transformation matrices, shape is (len(self), 3, 3).
        """
        if self.isscalar:
            return np.transpose(self.itrs2gcrs_dot)
        else:
            return np.transpose(self.itrs2gcrs_dot, axes=[0, 2, 1])

    @cache.property
    def jd_int(self):
        """Integer part of Julian Day

        Internally, Astropy-Time-objects are represented by two double precision numbers `jd1` and `jd2`. See [1,
        section Creating a Time Object]. Although, these typically will be the integer and fractional part of the
        Julian Day, there are no guarantees. An example is as follows:

            > import astropy.time
            > t1 = astropy.time.Time(val=2457617.5, val2=0., format='jd', scale='utc')
            > t1.jd1, t1.jd2
            (2457617.5, 0.0)
            > t2 = astropy.time.Time(val=t1.jd, val2=0., format='jd', scale='utc')
            > t2.jd1, t2.jd2
            (2457618.0, -0.5)

        To ensure consistency, we therefore add two properties `jd_int` and `jd_frac` where the integer part is
        guaranteed to be a "half-integer" (e.g. 2457617.5) and the fractional part is guaranteed to be a float in the
        range [0., 1.). The parts are calculated from `jd1` and `jd2` to preserve precision.

        Returns:
            Numpy-float scalar or array with (half-)integer part of Julian Day.
        """
        return self.jd1 - self._jd_delta

    @cache.property
    def jd_frac(self):
        """Fractional part of Julian Day

        See the docstring of `jd_int` for more information.

        Returns:
            Numpy-float scalar or array with fractional part of Julian Day, in the range [0., 1.).
        """
        return self.jd2 + self._jd_delta

    @cache.property
    def sec_of_day(self):
        """Seconds since midnight

        Note   -  Does not support leap seconds

        Returns:
            int/list: seconds since midnight

        """
        if self.isscalar:
            return self.datetime.hour * 60 * 60 + self.datetime.minute * 60 + self.datetime.second

        return np.array([d.hour * 60 * 60 + d.minute * 60 + d.second for d in self.datetime])

    @cache.property
    def yydddsssss(self):
        """Text representation YY:DDD:SSSSS

        YY     - decimal year without century
        DDD    - zero padded decimal day of year
        SSSSS  - zero padded seconds since midnight

        Note   -  Does not support leap seconds

        Returns:
            string/list: Time converted to yydddssss format

        """
        if self.isscalar:
            return self.datetime.strftime("%y:%j:") + str(self.sec_of_day).zfill(5)

        return [d.strftime("%y:%j:") + str(s).zfill(5) for d, s in zip(self.datetime, self.sec_of_day)]

    @cache.property
    def mean(self):
        """Mean time

        Returns:
            Time:    Time object containing the mean time
        """
        if self.isscalar:
            return self

        return Time(np.mean(self.utc.jd), scale="utc", format="jd")

    @cache.property
    def _jd_delta(self):
        """Delta between jd1 and jd_int

        This is a helper function used by `jd_int` and `jd_frac` to find the difference to `jd1` and `jd2`
        respectively. See the docstring of `jd_int` for more information.

        Returns:
            Numpy-float scalar or array with difference between `jd1` and the integer part of Julian Day.
        """
        return self.jd1 - (np.floor(self.jd - 0.5) + 0.5)

    @cache.property
    def mjd_int(self):
        """Integer part of Modified Julian Day

        In general, we have that MJD = JD - 2400000.5. See the docstring of `jd_int` for more information.

        Returns:
            Numpy-float scalar or array with the integer part of Modified Julian Day.
        """
        return self.jd_int - 2400000.5

    @cache.property
    def mjd_frac(self):
        """Fractional part of Modified Julian Day

        See the docstring of `jd_int` for more information. The way we have defined `jd_int` and `jd_frac` means that
        `mjd_frac` will be equal to `jd_frac`.

        Returns:
            Numpy-float scalar or array with the fractional part of Modified Julian Day, in the range [0., 1.).
        """
        return self.jd_frac

    @cache.property
    def gpsweek(self):
        """The GPS week for the given time-object

        A GPS week starts on Sunday. GPS week 0 starts on January 6, 1980 at 00:00 UTC, but follows the GPS time
        scale. That is, leap seconds are not included.

        Returns:
            Numpy-float scalar or array with GPS week.
        """
        gps_delta = self.gps - Time(val="1980-01-06 00:00:00", scale="utc", format="iso")
        return np.floor(gps_delta.jd / 7)

    @cache.property
    def gpsday(self):
        """The day of the week following GPS conventions (Sunday: 0, ..., Saturday: 6)

        The calculation is done by the datetime.weekday-method. However, this uses the convention Monday: 0, ...,
        Sunday: 6 so we correct for this by adding 1 modulo 7.

        Returns:
            Numpy-int scalar or array with the day of the week, Sunday: 0, ..., Saturday: 6.
        """
        if self.isscalar:
            return (self.gps.datetime.weekday() + 1) % 7
        else:
            gpsdays = np.empty(self.shape, dtype=int)
            gpsdays[:] = [(d.weekday() + 1) % 7 for d in self.gps.datetime]
            return gpsdays

    @cache.property
    def gpssec(self):
        """Number of seconds since the start of the GPS week

        A GPS week starts on Sunday at 00:00:00 GPS. GPS seconds are number of seconds since the start of the GPS
        week. This will be a number between 0 and 604800 (= 7 x 24 x 60 x 60).

        Returns:
            Numpy-float scalar or array with seconds since the start of the GPS week.
        """
        return (self.gps.gpsday + self.gps.jd_frac) * 86400

    @cache.function
    def sec_to_reference(self, ref_time):
        """Seconds since ref_time

        Args:
            ref_time (Time):  Reference epoch.

        Returns:
            float/array: Seconds since reference epoch
        """
        if not isinstance(ref_time, astropy.time.Time):
            ref_time = astropy.time.Time(ref_time, scale="utc")

        return (self - ref_time).sec

    @property
    def data(self):
        """Temporary warning about removing of data field

        Remove this method when all references to time.data are gone.
        """
        import sys
        from where.lib import log

        caller = sys._getframe(1)
        func_name = caller.f_code.co_name
        file_name = caller.f_code.co_filename
        line_num = caller.f_lineno
        log.dev("'time.data' is deprecated. Use 'time' instead in '{}' ({}:{})", func_name, file_name, line_num)

        return self

    def _get_delta_ut1_utc(self, jd1=None, jd2=None):
        """Calculate ut1 - utc based on EOP tables

        Reads the ut1 - utc field in the apriori EOP tables and interpolates for each time epoch.
        TODO: Take jd1, jd2 into account
        """
        from where import apriori

        eop = apriori.get("eop", time=self)

        return eop.ut1_utc

    def _get_delta_tdb_tt(self, jd1=None, jd2=None):
        """Calculate tdb - tt based on the IERS conventions

        The delta TDB - TT is calculated by way of TCG and TCB as described in the IERS Conventions, chapter 10.1
        [1]. See also the note `Calculation of Time Deltas in Where` available at `where/documents/notes/time_deltas`.

        """
        if hasattr(self, "_delta_tdb_tt"):
            return self._delta_tdb_tt

        if jd1 is None or jd2 is None:
            if self.scale in ("tt", "tdb"):
                jd1 = self._time.jd1
                jd2 = self._time.jd2
            else:
                raise ValueError("Accessing the delta_tdb_tt attribute is only possible for TT or TDB time scales")
        jd = jd1 + jd2

        # Conventional TCB - TCG (hf2002_iers wants jd in TDB, but TT ok)
        if self.isscalar:
            tcb_tcg = iers.hf2002_iers(jd)
        else:
            tcb_tcg = np.array([iers.hf2002_iers(t) for t in jd])

        if self.scale in ("tdb", "tcb"):  # jd is in TCB time scale
            # Equation (10.3) in [1]. Separate terms to avoid loss of precision
            tdb_tcb = (constant.TDB_0 - constant.L_B * (jd - constant.T_0) * 86400) / (1 - constant.L_B)

            # TCG - TT may be calculated following equation (10.1) in [1].
            tcg_tt = constant.L_G * (jd - constant.T_0) * 86400 - constant.L_G * (tdb_tcb + tcb_tcg)

        else:  # jd is in TT time scale
            # TCG - TT may be calculated following equation (10.1) in [1].
            tcg_tt = constant.L_G / (1 - constant.L_G) * (jd - constant.T_0) * 86400

            # Equation (10.3) in [1]. Separate terms to avoid loss of precision
            tdb_tcb = constant.TDB_0 - (constant.L_B * (tcb_tcg + tcg_tt) + constant.L_B * (jd - constant.T_0) * 86400)

        self._delta_tdb_tt = tdb_tcb + tcb_tcg + tcg_tt
        return self._delta_tdb_tt

    delta_tdb_tt = property(_get_delta_tdb_tt, astropy.time.Time._set_delta_tdb_tt)

    def __hash__(self):
        """Define a hash value

        We want arrays with the same time epochs to have equal hashes. This is done by applying a hash to a textual
        representation of the time array.
        """
        return hash(self.text_repr)

    def __eq__(self, other):
        """Compare if two Time objects are equal

        Use the text representation to say that two Time objects are equal if they have the exact same time
        epochs. (This was seemingly quite a bit faster than to do a `all(super().__eq__(other))`-call.)

        This was necessary to implement because things like Dict use __eq__ to double-check equality if two
        __hash__-values are the same (to guard against hash-collisions), and the __eq__ of Astropy-Time uses
        element-wise comparison (numpy-style).
        """
        if not isinstance(other, self.__class__):
            return False

        return self.text_repr == other.text_repr

    def __repr__(self):
        """Represent the time object with scale and format, and an indication of values
        """
        value = self.value if self.isscalar else "<{} epochs>".format(len(self))
        return (
            "{}(scale='{}', format='{}', value={})" "".format(self.__class__.__name__, self.scale, self.format, value)
        )


def _add_gps_scale():
    """Add scale with name gps

    We define the 'gps' time scale in terms of the 'tai' time scale. All conversions to other time scales are done by
    first converting to 'tai', and then converting from 'tai' to the given time scale. See the update of the MULTI_HOPS
    dictionary.

    The actual conversion to and from 'tai' is done by the local functions `_gpstai` and `_taigps`. We treat the 'gps'
    time scale as being at a constant offset of 19 seconds compared to 'tai'. However, GPS time is automatically
    steered to UTC on a daily basis. This means that our implementation of 'gps' is only accurate to a few hundred
    nano-seconds. See http://tycho.usno.navy.mil/gpstt.html for more details.

    Technical notes: This feels quite hacky and could easily break if astropy changes their implementation of time
    scales. Ideally, we should have been able to simply add the 'gps' scale in our subclass to the proper lists with
    some hooks for the conversion functions. Unfortunately, a lot of the information is kept in global tuples and dicts
    on the module level (core and formats), so we must hook everything in there instead. Most of the original
    information is in the first 100 lines of astropy/time/core.py in the Astropy source code. See [2].

    Furthermore, all conversion functions are assumed to be implemented by the ERFA library (essentially the SOFA
    library), so we need to monkey patch the astropy._erfa module as well :(

    There has been some discussion in the Astropy community on including 'gps' as a time scale (see for instance
    https://github.com/astropy/astropy/pull/1879), but for now the consensus seems to be that it will not be
    implemented. It might be useful for us to create a pull request with a time scale architecture that is easier to
    extend though.
    """
    # Rename the 'gps' format to 'gps_1980'
    astropy.time.core.TIME_FORMATS["gps_1980"] = astropy.time.core.TIME_FORMATS["gps"]
    del astropy.time.core.TIME_FORMATS["gps"]

    # Add 'gps' as a scale
    astropy.time.core.TIME_SCALES = ("gps",) + astropy.time.core.TIME_SCALES
    astropy.time.Time.SCALES = astropy.time.core.TIME_SCALES
    #    Time.SCALES = astropy.time.core.TIME_SCALES
    astropy.time.formats.TIME_SCALES = astropy.time.core.TIME_SCALES
    astropy.time.core.GEOCENTRIC_SCALES = ("gps",) + astropy.time.core.GEOCENTRIC_SCALES
    astropy.time.core.TIME_DELTA_TYPES = dict(
        (scale, scales)
        for scales in (
            astropy.time.core.GEOCENTRIC_SCALES,
            astropy.time.core.BARYCENTRIC_SCALES,
            astropy.time.core.ROTATIONAL_SCALES,
        )
        for scale in scales
    )
    astropy.time.core.TIME_DELTA_SCALES = sorted(astropy.time.core.TIME_DELTA_TYPES.keys())
    astropy.time.TimeDelta.SCALES = astropy.time.core.TIME_DELTA_SCALES
    astropy.time.formats.TIME_DELTA_SCALES = astropy.time.core.TIME_DELTA_SCALES

    # Add information about how to convert to and from 'gps'
    astropy.time.core.MULTI_HOPS.update(
        {
            ("gps", "tcb"): ("tai", "tt", "tdb"),
            ("gps", "tcg"): ("tai", "tt"),
            ("gps", "tdb"): ("tai", "tt"),
            ("gps", "tt"): ("tai",),
            ("gps", "ut1"): ("tai", "utc"),
            ("gps", "utc"): ("tai",),
        }
    )

    # Conversion functions between 'gps' and 'tai' are monkey patched into the astropy._erfa-module
    def _gpstai(gps1, gps2):
        """Convert from 'gps' time scale to 'tai'

        Args:
            gps1:  Float, part 1 of 'gps'-time as a two-part Julian Day.
            gps2:  Float, part 2 of 'gps'-time as a two-part Julian Day.

        Returns:
            2-tuple of floats representing 'tai' as a two-part Julian Day.
        """
        return gps1, gps2 + 19 / 86400

    def _taigps(tai1, tai2):
        """Convert from 'tai' time scale to 'gps'

        Args:
            gps1:  Float, part 1 of 'gps'-time as a two-part Julian Day.
            gps2:  Float, part 2 of 'gps'-time as a two-part Julian Day.

        Returns:
            2-tuple of floats representing 'tai' as a two-part Julian Day.
        """
        return tai1, tai2 - 19 / 86400

    astropy._erfa.gpstai = _gpstai
    astropy._erfa.taigps = _taigps


# Monkey patch to get 'gps' scale when this module is imported
_add_gps_scale()
