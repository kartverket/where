"""Try to detect clock breaks for a VLBI session

Description:
------------

Suspected clock breaks are added as events to the dataset.

"""
# Standard library imports

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import config
from where.lib import log
from where.data.time import Time


@plugins.register
def detect_clockbreaks(dset):
    """Try to detect a clock break

    The suspected clock breaks are added to the dataset as events. They will not be corrected automatically.

    TODO: Clean up code / better variable names, remove "magic" numbers, possibly put settings in config?
          Handle gaps in observations better, these are often misinterpreted as clock breaks
          These might be made faster by using a filter to test all observations at the same time?

    Args:
        dset (Dataset):  Information about model run.
    """
    log.info("Looking for clock breaks")

    clock_breaks = list()
    order_of_polynomial = config.tech.get("order_of_polynomial", section="vlbi_clock_poly", default=2).int

    for station in dset.unique("station"):
        # Merge together data where station is station 1 and station 2
        idx_1 = dset.filter(station_1=station)
        idx_2 = dset.filter(station_2=station)
        time = np.hstack((dset.time.utc.mjd[idx_1], dset.time.utc.mjd[idx_2])) - dset.time.utc.mjd[0]
        residual = np.hstack((dset.residual[idx_1], -dset.residual[idx_2]))

        # Make sure data are chronological
        idx_sort = np.argsort(time)
        time = time[idx_sort]
        residual = residual[idx_sort]

        # Add fields to dset for debug
        idx_site = np.hstack((np.where(idx_1)[0], np.where(idx_2)[0]))[idx_sort]
        dset.add_float(f"cb_{station}_residual", np.zeros(dset.num_obs), write_level="operational")
        dset[f"cb_{station}_residual"][idx_site] = residual
        dset.add_float(f"cb_{station}_value", np.zeros(dset.num_obs), write_level="detail")
        dset.add_float(f"cb_{station}_limit", np.zeros(dset.num_obs), write_level="detail")
        dset.add_float(f"cb_{station}_pred", np.zeros(dset.num_obs), write_level="detail")
        dset.add_float(f"cb_{station}_ratio", np.zeros(dset.num_obs), write_level="detail")

        # Test each observation for clock break
        start_obs = 0
        for obs in range(len(time) - 25):
            if obs - start_obs < 25:  # Need some observations to do polyfit
                continue

            # Fit a polynomial to the given data
            idx_fit = slice(np.maximum(start_obs, obs - 500), obs)  # Possibly better to limit on time instead of obs?
            p = np.polyfit(time[idx_fit], residual[idx_fit], order_of_polynomial)  # Same degree as clock correction
            poly = np.polyval(p, time[idx_fit])
            res = residual[idx_fit] - poly
            std_lim = 2 * np.std(res) * (1 + 4 * np.exp(-10 * ((obs - np.maximum(start_obs, obs - 500)) / 500) ** 2))
            # Gives higher limit when there are fewer observations in res

            # Test next observations
            model = np.polyval(
                p, time[obs + 1 : obs + 26]
            )  # Use many (=25) observations to avoid problems with outliers
            obs_res = residual[obs + 1 : obs + 26]
            dset[f"cb_{station}_value"][idx_site[obs]] = np.min(np.abs(obs_res - model))
            dset[f"cb_{station}_limit"][idx_site[obs]] = std_lim
            dset[f"cb_{station}_pred"][idx_site[obs]] = model[0]
            dset[f"cb_{station}_ratio"][idx_site[obs]] = np.min(np.abs(obs_res - model)) / std_lim

            # Register possible clock break
            if np.all(np.abs(obs_res - model) > std_lim):
                start_obs = np.min(np.where(time > time[obs])[0])  # Next epoch with observations
                time_cb = Time(dset.time.utc.mjd[0] + (time[obs] + time[start_obs]) / 2, fmt="mjd", scale="utc")
                clock_breaks.append((np.min(np.abs(obs_res - model)) / std_lim, time_cb, station))

    # Only actually add the biggest clock breaks, because big clock breaks creates smaller false clock breaks
    ratio_lim = max(cb[0] for cb in clock_breaks) / 3 if clock_breaks else 1
    for ratio, time, station in clock_breaks:
        if ratio > ratio_lim:
            dset.meta.add_event(time, "suspected_clock_break", station)
            cb_time = time.datetime.strftime(config.FMT_datetime)
            stars = "*" * int(np.ceil(np.log2(ratio)))
            log.check(f"Found possible clock break for {station} at {cb_time} ({stars})")
