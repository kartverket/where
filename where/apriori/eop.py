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
    > from where.data.time import Time
    > time = Time(np.arange(55000, 55002, 1e-3), fmt="mjd", scale="utc")

    > eop = apriori.get('eop', time=time)
    > eop.x
    array([ 0.09106713,  0.09107276,  0.0910784 , ...,  0.09697505,
            0.09697876,  0.09698251])



"""
# Standard library imports
from functools import lru_cache

# External library imports
from scipy import interpolate
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math import nputil

# Where imports
from where import parsers
from where.data.time import Time
from where.ext import iers_2010_wrapper as iers
from where.ext import hf_eop_wrapper as hf_eop
from where.lib import config
from where.lib import exceptions
from where.lib import log
from where.lib.unit import Unit

# Cache for EOP data read from file
_EOP_DATA = dict()


@plugins.register
def get_eop(time, models=None, pole_model=None, window=4, sources=None):
    """Get EOP data for the given time epochs

    Args:
        time (Time):   Time epochs for which to calculate EOPs.
        models (Tuple): Optional tuple of EOP models. If not given, the config setting is used.

    Returns:
        Eop: Object that calculates EOP corrections.

    """
    if not _EOP_DATA:

        sources = config.tech.get("eop_sources", value=sources).list
        for source in sources:
            _EOP_DATA.setdefault(source, {}).update(parsers.parse_key(file_key=f"eop_{source}").as_dict())

    return Eop(_EOP_DATA, time, models=models, pole_model=pole_model, window=window)


class Eop:
    """A class that can calculate EOP corrections.

    One instance of the `Eop`-class calculates corrections for a given set of time epochs (specified when the instance
    is created). However, all instances share a cache of results from the various functions calculating corrections.

    The class properties are also cached and since the data content of each Eop instance is static there is no need to
    reset the cache at any point.
    """

    _correction_cache = dict()

    def __init__(self, eop_data, time, models=None, pole_model=None, window=4, sources=None):
        """Create an Eop-instance that calculates EOP corrections for the given time epochs

        The interpolation window is based on https://hpiers.obspm.fr/iers/models/interp.f which uses 4 days.

        Args:
            eop_data (Dict): Dictionary of tabular EOP data, typically read from file.
            time (Time):     Time epochs for which to calculate EOPs.
            models (Tuple):  Optional tuple of EOP correction models. If not given, the config setting is used.
            window (Int):    Number of days to use as interpolation window.
        """
        if time.scale == "ut1":
            raise ValueError(f"Time scale of 'time' cannot be 'ut1'")
        self.window = window
        self.time = time
        self.sources = sources
        self.data = self.pick_data(eop_data, self.time, self.window, sources)
        self.calculate_leap_second_offset()

        # Figure out which correction models to use
        self.models = config.tech.eop_models.tuple if models is None else models

        if "rg_zont2" in self.models:
            self.remove_low_frequency_tides()

        # Figure out which pole model to use:
        self.pole_model = config.tech.get("eop_pole_model", value=pole_model, default=None).str
        if self.pole_model == "mean_2015":
            # Read the tabulated data needed for the model
            data = parsers.parse_key("eop_mean_pole_2015").as_dict()
            self.mean_pole_last_idx = len(data["year"]) - 1
            self.mean_pole_years = interpolate.interp1d(
                data["year"], data["year"], kind="previous", fill_value="extrapolate"
            )
            self.mean_pole_idx = interpolate.interp1d(
                data["year"], range(len(data["year"])), kind="previous", fill_value=np.nan, bounds_error=False
            )
            self.mean_pole_x = interpolate.interp1d(
                range(len(data["x"])), data["x"], kind="previous", fill_value="extrapolate"
            )
            self.mean_pole_y = interpolate.interp1d(
                range(len(data["y"])), data["y"], kind="previous", fill_value="extrapolate"
            )

    @staticmethod
    def pick_data(eop_data, time, window, sources):
        """Pick out subset of eop_data relevant for the given time epochs and interpolation window

        Args:
            eop_data (Dict):   Dictionary of EOP data indexed by MJD dates.
            time (Time):       Time epochs for which to calculate EOPs.
            window (Int):      Interpolation window [days].

        Returns:
            Dict: EOP data subset to the time period needed.
        """
        if time.size == 1:
            start_time = np.floor(time.utc.mjd) - window // 2
            end_time = np.ceil(time.utc.mjd) + window // 2
        else:
            start_time = np.floor(time.utc.mjd.min()) - window // 2
            end_time = np.ceil(time.utc.mjd.max()) + window // 2

        sources = config.tech.get("eop_sources", value=sources).list
        for source in sources:
            try:
                picked_data = {d: eop_data[source][d].copy() for d in np.arange(start_time, end_time + 1)}
                eop_path = config.files.path(f"eop_{source}")
                log.debug(f"Using a priori EOP values from {eop_path} ")
                return picked_data
            except KeyError:
                pass

        # No data found if we reached this point
        paths = [str(config.files.path(f"eop_{k}")) for k in sources]
        raise exceptions.MissingDataError(
            "Not all days in the time period {:.0f} - {:.0f} MJD were found in EOP-files {}"
            "".format(start_time, end_time, ", ".join(paths))
        )

    def calculate_leap_second_offset(self):
        """Calculate leap second offsets for each day

        Use the difference between UTC and TAI as a proxy for the leap second offset. The leap second offset is
        calculated and stored to the EOP data-dictionary. This is used to correct for the leap second jumps when
        interpolating the UT1 - UTC values.
        """
        days = Time(np.array(list(self.data.keys())), fmt="mjd", scale="utc")
        leap_offset = np.round((days.utc.mjd - days.tai.mjd) * Unit.day2seconds)
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
            t = Time(mjd, fmt="mjd", scale="utc")
            self.data[mjd]["ut1_utc"] -= iers.rg_zont2(t)[0]

    @property
    @lru_cache()
    @Unit.register("arcseconds")
    def x(self):
        """X-motion of the Celestial Intermediate Pole

        See section 5.5.1 in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("x")
        values += self._corrections(
            ("ortho_eop", iers.ortho_eop, 0, 1e-6),
            ("pmsdnut2", iers.pmsdnut2, 0, 1e-6),
            ("hf_eop_xyu", hf_eop.hf_eop_xyu, 0, 1e-6),
        )
        return values

    @property
    @lru_cache()
    @Unit.register("arcseconds per day")
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

    @property
    @lru_cache()
    @Unit.register("arcseconds")
    def x_pole(self):
        return getattr(self, f"x_{self.pole_model}")()

    #
    # x-pole models
    #
    @Unit.register("arcseconds")
    def x_secular(self):
        """Returns the x-coordinate of the secular pole

        See chapter 7 in IERS Conventions, ::cite:`iers2010`:

        Returns:
            Array: Secular X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        # IERS conventions 2010 v.1.2.0 (chapter 7, equation 21)
        return (55.0 + 1.677 * (self.time.jyear - 2000)) * Unit.milliarcsec2arcsec

    @Unit.register("arcseconds")
    def x_mean_2015(self):
        """x-coordindate of Conventional mean pole model version 2015

        Reimplementation of IERS Conventions 2010 Software function IERS_CMP_2015.F
        (ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/software/IERS_CMP_2015.F)

        Units: Arcseconds
        """
        epochs = self.time.jyear
        mean_pole_idx = self.mean_pole_idx(epochs)
        # Inside of tabulated data range
        in_range_idx = ~np.isnan(mean_pole_idx)
        dt = epochs[in_range_idx] - self.mean_pole_years(epochs[in_range_idx])
        idx = mean_pole_idx[in_range_idx]
        x = np.full(len(epochs), fill_value=np.nan)
        x[in_range_idx] = self.mean_pole_x(idx) + dt * (self.mean_pole_x(idx + 1) - self.mean_pole_x(idx))

        # Extrapolate outside of tabulated data range
        dt = epochs[~in_range_idx] - self.mean_pole_years(epochs[~in_range_idx])
        x[~in_range_idx] = self.mean_pole_x(self.mean_pole_last_idx) + dt * (
            self.mean_pole_x(self.mean_pole_last_idx) - self.mean_pole_x(self.mean_pole_last_idx - 1)
        )
        return x

    @Unit.register("arcseconds")
    def x_mean_2010(self):
        """x-coordindate of Conventional mean pole model version 2010

        Reimplementation of IERS Conventions 2010 Software function IERS_CMP_2015.F
        (ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/software/IERS_CMP_2015.F)

        Units: Arcseconds
        """
        epochs = self.time.jyear
        dt = epochs - 2000.0
        idx = dt < 10
        x = np.zeros(len(epochs))
        x[idx] = 0.055_974 + 0.001_824_3 * dt[idx] + 0.000_184_13 * dt[idx] ** 2 + 0.000_007_024 * dt[idx] ** 3
        x[~idx] = 0.23513 + 0.007_614_1 * dt[~idx]
        return x

    @Unit.register("arcseconds")
    def x_mean_2003(self):
        """x-coordindate of Conventional mean pole model version 2003

        Reimplementation of IERS Conventions 2010 Software function IERS_CMP_2015.F
        (ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/software/IERS_CMP_2015.F)

        Units: Arcseconds
        """
        return 0.054 + 0.00083 * (self.time.jyear - 2000.0)

    @property
    @lru_cache()
    @Unit.register("arcseconds")
    def y_pole(self):
        return getattr(self, f"y_{self.pole_model}")()

    #
    # y-pole models
    #
    @Unit.register("arcseconds")
    def y_secular(self):
        """Returns the x-coordinate of the secular pole

        See chapter 7 in IERS Conventions, ::cite:`iers2010`:

        Returns:
            Array: Mean X-motion of the CIP, one value for each time epoch [arcseconds].
        """
        # IERS conventions 2010 v.1.2.0 (chapter 7, equation 21)
        return (320.5 + 3.460 * (self.time.jyear - 2000)) * Unit.milliarcsec2arcsec

    @Unit.register("arcseconds")
    def y_mean_2015(self):
        """y-coordindate of Conventional mean pole model version 2015

        Reimplementation of IERS Conventions 2010 Software function IERS_CMP_2015.F
        (ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/software/IERS_CMP_2015.F)

        Units: Arcseconds
        """
        epochs = self.time.jyear
        mean_pole_idx = self.mean_pole_idx(epochs)
        # Inside of tabulated data range
        in_range_idx = ~np.isnan(mean_pole_idx)
        dt = epochs[in_range_idx] - self.mean_pole_years(epochs[in_range_idx])
        idx = mean_pole_idx[in_range_idx]
        y = np.full(len(epochs), fill_value=np.nan)
        y[in_range_idx] = self.mean_pole_y(idx) + dt * (self.mean_pole_y(idx + 1) - self.mean_pole_y(idx))

        # Extrapolate outside of tabulated data range
        dt = epochs[~in_range_idx] - self.mean_pole_years(epochs[~in_range_idx])
        y[~in_range_idx] = self.mean_pole_y(self.mean_pole_last_idx) + dt * (
            self.mean_pole_y(self.mean_pole_last_idx) - self.mean_pole_y(self.mean_pole_last_idx - 1)
        )
        return y

    @Unit.register("arcseconds")
    def y_mean_2010(self):
        """y-coordindate of Conventional mean pole model version 2010

        Reimplementation of IERS Conventions 2010 Software function IERS_CMP_2015.F
        (ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/software/IERS_CMP_2015.F)

        Units: Arcseconds
        """
        epochs = self.time.jyear
        dt = epochs - 2000.0
        idx = dt < 10
        y = np.zeros(len(epochs))
        y[idx] = 0.346_346 + 0.001_789_6 * dt[idx] - 0.000_107_29 * dt[idx] ** 2 - 0.000_000_908 * dt[idx] ** 3
        y[~idx] = 0.358_891 - 0.000_628_7 * dt[~idx]
        return y

    @Unit.register("arcseconds")
    def y_mean_2003(self):
        """y-coordindate of Conventional mean pole model version 2003

        Reimplementation of IERS Conventions 2010 Software function IERS_CMP_2015.F
        (ftp://maia.usno.navy.mil/conventions/2010/2010_update/chapter7/software/IERS_CMP_2015.F)

        Units: Arcseconds
        """
        return 0.357 + 0.00395 * (self.time.jyear - 2000.0)

    @property
    @lru_cache()
    @Unit.register("arcseconds")
    def y(self):
        """Y-motion of the Celestial Intermediate Pole

        See section 5.5.1 in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: Y-motion of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("y")
        values += self._corrections(
            ("ortho_eop", iers.ortho_eop, 1, 1e-6),
            ("pmsdnut2", iers.pmsdnut2, 1, 1e-6),
            ("hf_eop_xyu", hf_eop.hf_eop_xyu, 1, 1e-6),
        )
        return values

    @property
    @lru_cache()
    @Unit.register("arcseconds per day")
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

    @property
    @lru_cache()
    @Unit.register("seconds")
    def ut1_utc(self):
        """Delta between UT1 and UTC

        See section 5.5.3 in IERS Conventions, :cite:`iers2010`. Does correction for leap second jumps before
        interpolation.

        Reapplies low frequency tides if these were removed before interpolation.

        Returns:
            Array: UT1 - UTC, one value for each time epoch [seconds].
        """
        values = self._interpolate_table("ut1_utc", leap_second_correction=True)
        values += self._corrections(
            ("ortho_eop", iers.ortho_eop, 2, 1e-6),
            ("utlibr", iers.utlibr, 0, 1e-6),
            ("hf_eop_xyu", hf_eop.hf_eop_xyu, 2, 1e-6),
            ("rg_zont2", iers.rg_zont2, 0, 1),
        )

        # low frequency tides
        # if "rg_zont2" in self.models:
        #    values += nputil.take(iers.rg_zont2(self.time), 0) # Column 0 contains ut1_utc corrections
        return values

    @property
    @lru_cache()
    @Unit.register("seconds per day")
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

        return values

    @property
    @lru_cache()
    @Unit.register("seconds")
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

    @property
    @lru_cache()
    @Unit.register("arcseconds")
    def dx(self):
        """X-offset of the Celestial Intermediate Pole

        See section 5.5.? in IERS Conventions, :cite:`iers2010`.

        Returns:
            Array: X-offset of the CIP, one value for each time epoch [arcseconds].
        """
        values = self._interpolate_table("dx")
        return values

    @property
    @lru_cache()
    @Unit.register("arcseconds")
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
        offsets = range(int(-np.ceil(self.window / 2) + 1), int(np.floor(self.window / 2) + 1))

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

        if self.time.size == 1:
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
        corrections = 0 if self.time.size == 1 else np.zeros(self.time.size)
        for name, correction_func, out_idx, factor in correction_models:
            if name not in self.models:
                continue

            corrections += factor * nputil.take(correction_func(self.time), out_idx)
        return corrections

    # Add methods to deal with units for Eop-properties (set by @Unit.register)
    convert_to = Unit.convert_factory(__name__)
    unit_factor = staticmethod(Unit.factor_factory(__name__))
    unit = staticmethod(Unit.unit_factory(__name__))
