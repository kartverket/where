"""Get apriori data for earth orientation parameters

Description:
------------

Reads data for earth orientation parameters from the two data files (eopc04_IAU2000.62-now and eopc04_extended.dat). If
there are data in both files, the former is prioritized.

The following parameters are provided:

========== ============================================================
Parameter  Description
========== ============================================================
x          X-polar motion pole in arcseconds
y          Y-polar motion pole in arcseconds
ut1_utc    UT1 - UTC in seconds
lod        Length of day in seconds
dx         X-celestial pole offset in arcseconds
dy         Y-celestial pole offset in arcseconds
========== ============================================================

The tabular values are interpolated and high frequency tides are added to the interpolated values based on section 5.5
in the IERS Conventions, :cite:`iers2010`.

Equation (5.2) in the IERS Conventions indicates that the appropriate time scale for the corrections is Terrestrial
Time (TT).

Example:
--------

To use the earth orientation parameters, simply get them from the apriori-package specifying the time epochs you need::

    eop = apriori.get('eop', time=time)

After the Eop-object is created, you can get the different values as properties::

    eop.x
    eop.y
    eop.ut1_utc
    eop.lod
    eop.dx
    eop.dy

Each of these returns a numpy array with one value for each time epoch.  A fully working example might look as
follows::

    > import numpy as np
    > from where import apriori
    > from where.lib.time import Time
    > time = Time(np.arange(55000, 55002, 1e-3), format='mjd')

    > eop = apriori.get('eop', time=time)
    > eop.x
    array([ 0.09106713,  0.09107276,  0.0910784 , ...,  0.09697505,
            0.09697876,  0.09698251])



"""
# Standard library imports
import math

# External library imports
from scipy import interpolate
import numpy as np

# Where imports
from where.lib import cache
from where.lib import config
from where.lib import files
from where.ext import iers_2010 as iers
from where import parsers
from where.lib import plugins
from where.lib.exceptions import MissingDataError
from where.lib.time import Time
from where.lib.unit import unit

# Cache for EOP data read from file
_EOP_DATA = dict()

# List of files to read for EOP data, last file is prioritized for overlapping data
_EOP_FILE_KEYS = {"c04": ("eop_c04_extended", "eop_c04"), "bulletin_a": ("eop_bulletin_a",)}


@plugins.register
def get_eop(time, models=None, window=4, source=None):
    """Get EOP data for the given time epochs

    Read EOP data from the eopc04-files. Both files are read, with data from the regular IAU file
    (eopc04_IAU2000.62-now) being prioritized, see `_EOP_FILE_KEYS`.

    Args:
        time (Time):   Time epochs for which to calculate EOPs.
        models (Tuple): Optional tuple of EOP models. If not given, the config setting is used.

    Returns:
        Eop: Object that calculates EOP corrections.

    """
    # Read the extended and the regular EOP data file (overlapping dates are overwritten by the latter)
    if not _EOP_DATA:
        source = config.tech.get("eop_source", value=source).str
        for file_key in _EOP_FILE_KEYS[source]:
            _EOP_DATA.update(parsers.parse_key(file_key=file_key).as_dict())

    return Eop(_EOP_DATA, time, models=models, window=window)


class Eop:
    """A class that can calculate EOP corrections.

    One instance of the `Eop`-class calculates corrections for a given set of time epochs (specified when the instance
    is created). However, all instances share a cache of results from the various functions calculating corrections.
    """

    _correction_cache = dict()

    def __init__(self, eop_data, time, models=None, window=4):
        """Create an Eop-instance that calculates EOP corrections for the given time epochs

        The interpolation window is based on https://hpiers.obspm.fr/iers/models/interp.f which uses 4 days.

        Args:
            eop_data (Dict): Dictionary of tabular EOP data, typically read from file.
            time (Time):     Time epochs for which to calculate EOPs.
            models (Tuple):  Optional tuple of EOP correction models. If not given, the config setting is used.
            window (Int):    Number of days to use as interpolation window.
        """
        self.window = window
        self.time = time
        self.data = self.pick_data(eop_data, self.time, self.window)
        self.calculate_leap_second_offset()

        # Figure out which correction models to use
        self.models = config.tech.get("eop_models", value=models, default="").list

        if "rg_zont2" in self.models:
            self.remove_low_frequency_tides()

        self._mean_pole_cache = dict()

    @staticmethod
    def pick_data(eop_data, time, window):
        """Pick out subset of eop_data relevant for the given time epochs and interpolation window

        Args:
            eop_data (Dict):   Dictionary of EOP data indexed by MJD dates.
            time (Time):       Time epochs for which to calculate EOPs.
            window (Int):      Interpolation window [days].

        Returns:
            Dict: EOP data subset to the time period needed.
        """
        if time.isscalar:
            start_time = np.floor(time.utc.mjd) - window // 2
            end_time = np.ceil(time.utc.mjd) + window // 2
        else:
            start_time = np.floor(time.min().utc.mjd) - window // 2
            end_time = np.ceil(time.max().utc.mjd) + window // 2

        try:
            return {d: eop_data[d].copy() for d in np.arange(start_time, end_time + 1)}
        except KeyError:
            paths = [str(files.path(k)) for k in _EOP_FILE_KEYS]
            raise MissingDataError(
                "Not all days in the time period {:.0f} - {:.0f} MJD were found in EOP-files {}"
                "".format(start_time, end_time, ", ".join(paths))
            )

    def calculate_leap_second_offset(self):
        """Calculate leap second offsets for each day

        Use the difference between UTC and TAI as a proxy for the leap second offset. The leap second offset is
        calculated and stored to the EOP data-dictionary. This is used to correct for the leap second jumps when
        interpolating the UT1 - UTC values.
        """
        days = Time(np.array(list(self.data.keys())), format="mjd", scale="utc")
        leap_offset = np.round((days.utc.mjd - days.tai.mjd) * unit.day2seconds)
        daily_offset = {int(d): lo for d, lo in zip(days.mjd, leap_offset)}

        for d, lo in daily_offset.items():
            self.data[d]["leap_offset"] = lo

    def remove_low_frequency_tides(self):
        """Remove the effect of low frequency tides.

        Tidal variations in the Earth's rotation with periods from 5 days to 18.6 years is present in the UT1-UTC time
        series as described in the IERS Conventions 2010 chapter 8.1. To improve the interpolation of the UT1-UTC time
        series this effect can be removed. In that case the effect needs to be added again to the final interpolated
        values.
        """
        for mjd in self.data.keys():
            # Julian centuries since J2000
            t = Time(mjd, format="mjd")
            t_julian_centuries = (t.tt.jd - 2451545.0) / 36525
            dut1_corr = iers.rg_zont2(t_julian_centuries)[0]
            self.data[mjd]["ut1_utc"] -= dut1_corr

    @cache.property
    @unit.register("arcseconds")
    def x(self):
        """X-motion of the Celestial Intermediate Pole

        See section 5.5.1 in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("x")
        values += self._corrections(("ortho_eop", iers.ortho_eop, 0, 1e-6), ("pmsdnut2", iers.pmsdnut2, 0, 1e-6))
        return values

    @cache.property
    @unit.register("arcseconds per day")
    def x_rate(self):
        """X-motion of the Celestial Intermediate Pole

        See section 5.5.1 in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("x", derivative_order=1)
        # values += self._corrections(('ortho_eop', iers.ortho_eop, 0, 1e-6),
        #                            ('pmsdnut2', iers.pmsdnut2, 0, 1e-6))
        return values

    @cache.property
    @unit.register("arcseconds")
    def x_mean(self):
        """Returns the the mean value for the X coordinate of the Celestial Intermediate Pole

        See chapter 7 in IERS Conventions, ::cite:`iers2010`:

        Returns:
            Array: Mean X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        return self._get_mean_pole("x")

    @cache.property
    @unit.register("arcseconds")
    def y(self):
        """Y-motion of the Celestial Intermediate Pole

        See section 5.5.1 in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: Y-motion of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("y")
        values += self._corrections(("ortho_eop", iers.ortho_eop, 1, 1e-6), ("pmsdnut2", iers.pmsdnut2, 1, 1e-6))
        return values

    @cache.property
    @unit.register("arcseconds per day")
    def y_rate(self):
        """X-motion of the Celestial Intermediate Pole

        See section 5.5.1 in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("y", derivative_order=1)
        # values += self._corrections(('ortho_eop', iers.ortho_eop, 0, 1e-6),
        #                            ('pmsdnut2', iers.pmsdnut2, 0, 1e-6))
        return values

    @cache.property
    @unit.register("arcseconds")
    def y_mean(self):
        """Returns the the mean value for the X coordinate of the Celestial Intermediate Pole

        See chapter 7 in IERS Conventions, ::cite:`iers2010`:

        Returns:
            Array: Mean X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        return self._get_mean_pole("y")

    @cache.property
    @unit.register("seconds")
    def ut1_utc(self):
        """Delta between UT1 and UTC

        See section 5.5.3 in IERS Conventions, :cite:`iers2010`. Does correction for leap second jumps before
        interpolation.

        Reapplies low frequency tides if these were removed before interpolation.

        Returns:
            Array: UT1 - UTC, one value for each time epoch [seconds].
        """
        values = self._interpolate_table("ut1_utc", leap_second_correction=True)
        values += self._corrections(("ortho_eop", iers.ortho_eop, 2, 1e-6), ("utlibr", iers.utlibr, 0, 1e-6))

        # low frequency tides
        if "rg_zont2" in self.models:
            correction_cache = self._correction_cache.setdefault("rg_zont2", dict())
            # Julian centuries since J2000
            t_julian_centuries = (self.time.tt.jd - 2451545.0) / 36525

            if self.time.isscalar:
                mjd = self.time.tt.mjd
                if mjd not in correction_cache:
                    correction_cache[mjd] = iers.rg_zont2(t_julian_centuries)[0]
                dut1_corr = correction_cache[mjd]
            else:
                dut1_corr = list()
                for t in self.time.tt:
                    if t.mjd not in correction_cache:
                        t_julian_centuries = (t.tt.jd - 2451545.0) / 36525
                        correction_cache[t.mjd] = iers.rg_zont2(t_julian_centuries)[0]
                    dut1_corr.append(correction_cache[t.mjd])

            values += dut1_corr
        return values

    @cache.property
    @unit.register("seconds per day")
    def ut1_utc_rate(self):
        """Delta between UT1 and UTC

        See section 5.5.3 in IERS Conventions, :cite:`iers2010`. Does correction for leap second jumps before
        interpolation.

        Reapplies low frequency tides if these were removed before interpolation.

        TODO: apply models based on eop.models
        Only works if eop.models = ()

        Returns:
            Array: UT1 - UTC, one value for each time epoch [seconds].
        """
        values = self._interpolate_table("ut1_utc", leap_second_correction=True, derivative_order=1)
        # values += self._corrections(("ortho_eop", iers.ortho_eop, 2, 1e-6), ("utlibr", iers.utlibr, 0, 1e-6))

        # Low frequency tides
        #         if "rg_zont2" in self.models:
        #             correction_cache = self._correction_cache.setdefault("rg_zont2", dict())
        #             # Julian centuries since J2000
        #             t_julian_centuries = (self.time.tt.jd - 2451545.0) / 36525
        #
        #             if self.time.isscalar:
        #                 mjd = self.time.tt.mjd
        #                 if mjd not in correction_cache:
        #                     correction_cache[mjd] = iers.rg_zont2(t_julian_centuries)[0]
        #                 dut1_corr = correction_cache[mjd]
        #             else:
        #                 dut1_corr = list()
        #                 for t in self.time.tt:
        #                     if t.mjd not in correction_cache:
        #                         t_julian_centuries = (t.tt.jd - 2451545.0) / 36525
        #                         correction_cache[t.mjd] = iers.rg_zont2(t_julian_centuries)[0]
        #                     dut1_corr.append(correction_cache[t.mjd])
        #
        #             values += dut1_corr
        #         return values
        return values

    @cache.property
    @unit.register("seconds")
    def lod(self):
        """Length of day

        See section 5.5.3 in IERS Conventions, :cite:`iers2010`.

        TODO: How should this be implemented? Previous implementation simply used nearest value. Is this in the
        conventions?

        Returns:
            Array: Length of day, one value for each time epoch [seconds].
        """
        values = self._interpolate_table("lod")
        return values

    @cache.property
    @unit.register("arcseconds")
    def dx(self):
        """X-offset of the Celestial Intermediate Pole

        See section 5.5.? in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: X-offset of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("dx")
        return values

    @cache.property
    @unit.register("arcseconds")
    def dy(self):
        """Y-offset of the Celestial Intermediate Pole

        See section 5.5.? in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: Y-offset of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("dy")

        return values

    def _interpolate_table(self, key, leap_second_correction=False, derivative_order=0):
        """Interpolate daily values to the given time epochs

        Uses Lagrange interpolation with the given interpolation window.

        We have observed that the Lagrange interpolation introduces instabilities when the EOP data are constant (as
        for instance in the VASCC-data). In this case, we force the Lagrange polynomial to be constant.

        Args:
            key (String):                   Name of data to be interpolated, key in `self.data`.
            leap_second_correction (Bool):  Whether data should be corrected for leap seconds before interpolation.

        Returns:
            Array: Interpolated values, one value for each time epoch.
        """
        days = np.unique(self.time.utc.mjd_int)
        offsets = range(-math.ceil(self.window / 2) + 1, math.floor(self.window / 2) + 1)

        if leap_second_correction:
            leap = {
                d: np.array(
                    [self.data[d + o].get("leap_offset", np.nan) - self.data[d]["leap_offset"] for o in offsets]
                )
                for d in days
            }
            for lo in leap.values():
                lo[np.isnan(lo)] = 0
        else:
            leap = {d: 0 for d in days}

        table_values = {d: np.array([self.data[d + o][key] for o in offsets]) + leap[d] for d in days}
        interpolators = {d: interpolate.lagrange(offsets, v) for d, v in table_values.items()}
        for poly in interpolators.values():
            poly.c[np.abs(poly.c) < 1e-15] = 0  # Avoid numerical instabilities for constant values

        if derivative_order:
            interp_values = {
                d: np.polyder(ip, derivative_order)(self.time.utc.mjd_frac) for d, ip in interpolators.items()
            }
        else:
            interp_values = {d: ip(self.time.utc.mjd_frac) for d, ip in interpolators.items()}

        if self.time.isscalar:
            return interp_values[self.time.utc.mjd_int]

        values = np.empty(self.time.size)
        for day in days:
            idx = self.time.utc.mjd_int == day
            values[idx] = interp_values[day][idx]

        return values

    def _corrections(self, *correction_models):
        """Calculate corrections to tabular values

        The correction models are specified as tuples with name, function, output column and scale factor. Calls to the
        correction functions are cached since some correction functions are used by several EOP-values.

        Args:
            correction_models (Tuple): Specification of correction models (see above)

        Returns:
            Array: Corrections to tabular values, one value for each time epoch.
        """
        corrections = 0 if self.time.isscalar else np.zeros(self.time.size)
        for name, correction_func, out_idx, factor in correction_models:
            if name not in self.models:
                continue

            correction_cache = self._correction_cache.setdefault(name, dict())
            if self.time.isscalar:
                mjd = self.time.tt.mjd
                if mjd not in correction_cache:
                    correction_cache[mjd] = correction_func(mjd)
                corrections += factor * correction_cache[mjd][out_idx]
            else:
                for idx, mjd in enumerate(self.time.tt.mjd):
                    if mjd not in correction_cache:
                        correction_cache[mjd] = correction_func(mjd)

                    corrections[idx] += factor * correction_cache[mjd][out_idx]

        return corrections

    def _get_mean_pole(self, coord):
        """Calculate mean pole and cache results

        coord:        'x' or 'y'

        Returns:
            Array:     mean pole values [arcseconds]
        """
        version = config.tech.mean_pole_version.str
        key = coord + "_" + str(version)
        if key not in self._mean_pole_cache:
            mean_xp = np.empty(self.time.size)
            mean_yp = np.empty(self.time.size)
            # Calculate correction
            for obs, time in enumerate(self.time.tt):
                # Equation (7.25) IERS Conventions 2010
                mean_xp[obs], mean_yp[obs], _ = iers.iers_cmp_2015(version, time.jyear)
            self._mean_pole_cache["x_" + str(version)] = mean_xp
            self._mean_pole_cache["y_" + str(version)] = mean_yp
        return self._mean_pole_cache[key]

    # Add methods to deal with units for Eop-properties (set by @unit.register)
    convert_to = unit.convert_factory(__name__)
    unit_factor = staticmethod(unit.factor_factory(__name__))
    unit = staticmethod(unit.unit_factory(__name__))
