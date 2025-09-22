"""Calculate the delay caused by the thermal deformation

Description:
------------

Calculate the delay caused by the thermal deformation.

References:
-----------

.. [1] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
       IERS Technical Note No. 36, BKG (2010).
       http://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36.html

.. [2] http://hpiers.obspm.fr/combinaison/documentation/articles/Thermal_Expansion_Modelling_Radio_Telescopes_Nothnagel.pdf



"""
# External library imports
import numpy as np
from scipy import optimize

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import log


@plugins.register
def vlbi_thermal_deformation(dset):
    r"""Calculate thermal deformation at both stations

    According to Nothnagel [2]_ the calculated quantity measures how much the antenna axis offset brings the receiver
    closer to the incoming wavefront. Thus, the number is the opposite of a delay. This should be handled by
    calculating :math:`\delta \tau_1 - \delta \tau_2` instead of :math:`\delta \tau_2 - \delta \tau_1` as usual. We do
    this by subtracting the ``multiplier`` terms instead of adding.

    Args:
        dset (Dataset): Model data.

    Returns:
        Numpy array: Corrections in meters for each observation

    """
    temp_funcs = calculate_temperature_functions(dset)

    data_out = np.zeros(dset.num_obs)
    for multiplier in dset.for_each_suffix("station"):
        # Note the sign
        data_out -= multiplier * thermal_deformation_station(dset, temp_funcs)

    return data_out


def thermal_deformation_station(dset, temp_funcs):
    """Calculate thermal deformation at one station

    The foundation is assumed to be made of concrete and the antenna is assumed to be made of steel.

    Args:
        dset:          A Dataset containing model data.
        antenna_info:  Antenna info from Nothnagel's antenna info file.
        temp:          Array of temperature sine functions

    Returns:
        Numpy array with delay caused by thermal deformation in meters.
    """
    antenna_info = apriori.get("vlbi_antenna_info")

    # time delay for antenna and foundation
    dt_a = 2 / 24  # unit.hours2days
    dt_f = 6 / 24  # unit.hours2days
    delays = np.zeros(dset.num_obs)

    sin_a = np.sin(dset.site_pos.azimuth)
    cos_a = np.cos(dset.site_pos.azimuth)
    sin_e = np.sin(dset.site_pos.elevation)
    cos_e = np.cos(dset.site_pos.elevation)
    cos_d = np.cos(dset.src_dir.declination)

    for ivsname in dset.unique("ivsname"):
        site_id = dset.meta["station"][ivsname]["site_id"] if ivsname in dset.meta else ""
        if ivsname not in antenna_info:
            log.warn(f"Missing antenna specifications for ivsname {str(ivsname)!r} ({str(site_id)}). Thermal deformation correction set to zero")
            continue

        idx = dset.filter(ivsname=ivsname)
        AO = antenna_info[ivsname]["axis_offset"]
        axis_type = antenna_info[ivsname]["mount"]
        focus_type = antenna_info[ivsname]["focus"]
        gamma_f = antenna_info[ivsname]["coefficient_foundation"]
        gamma_a = antenna_info[ivsname]["coefficient_fixed_axis"]
        T_0 = antenna_info[ivsname]["reference_temperature"]
        h_f = antenna_info[ivsname]["height_foundation"]
        h_p = antenna_info[ivsname]["fixed_axis"]
        h_v = antenna_info[ivsname]["distance_antenna_vertex"]
        h_s = antenna_info[ivsname]["height_focus"]
        T = temp_funcs[ivsname]
        t = dset.time.utc.mjd

        if focus_type == "FO_PRIM":
            F_a = 0.9
        elif focus_type == "FO_SECN":
            F_a = 1.8
        else:
            log.warn(f"Unknown antenna focus type {focus_type!r} for {ivsname} ({site_id}). Correction set to zero")
            continue

        if axis_type == "MO_AZEL":
            delays[idx] = gamma_f * (T(t - dt_f)[idx] - T_0) * (h_f * sin_e[idx]) + gamma_a * (
                T(t - dt_a)[idx] - T_0
            ) * (h_p * sin_e[idx] + AO * cos_e[idx] + h_v - F_a * h_s)
        elif axis_type == "MO_EQUA":
            delays[idx] = gamma_f * (T(t - dt_f)[idx] - T_0) * (h_f * sin_e[idx]) + gamma_a * (
                T(t - dt_a)[idx] - T_0
            ) * (h_p * sin_e[idx] + AO * cos_d[idx] + h_v - F_a * h_s)
        elif axis_type == "MO_XYNO":
            delays[idx] = gamma_f * (T(t - dt_f)[idx] - T_0) * (h_f * sin_e[idx]) + gamma_a * (
                T(t - dt_a)[idx] - T_0
            ) * (h_p * sin_e[idx] + AO * np.sqrt(1 - (cos_e[idx] * cos_a[idx]) ** 2) + h_v - F_a * h_s)
        elif axis_type == "MO_XYEA":
            delays[idx] = gamma_f * (T(t - dt_f)[idx] - T_0) * (h_f * sin_e[idx]) + gamma_a * (
                T(t - dt_a)[idx] - T_0
            ) * (h_p * sin_e[idx] + AO * np.sqrt(1 - (cos_e[idx] * sin_a[idx]) ** 2) + h_v - F_a * h_s)
        elif axis_type == "MO_RICH":
            # Special case for RICHMOND station
            phi = np.radians(39.06)
            dlam = np.radians(-0.12)
            delays[idx] = gamma_f * (T(t - dt_f)[idx] - T_0) * (h_f * sin_e[idx]) + gamma_a * (
                T(t - dt_a)[idx] - T_0
            ) * (
                h_p * sin_e[idx]
                + AO
                * np.sqrt(
                    1
                    - (
                        sin_e[idx] * np.sin(phi)
                        + cos_e[idx] * np.cos(phi) * (cos_a[idx] * np.cos(dlam) + sin_a[idx] * np.sin(dlam))
                    )
                    ** 2
                )
                + h_v
                - F_a * h_s
            )
        else:
            log.warn(f"Unknown antenna axis type {axis_type!r} for {ivsname} ({site_id}). Correction set to zero")
            continue

        if any(np.isnan(delays)):
            delays[idx] = 0
    return delays


def calculate_temperature_functions(dset):
    """Estimates a temperature sine function based on the air temperature recorded at each site

    If the number of data points is too small to estimate a sine function a simple constant function based on the mean
    value of the datapoints is used instead.

    Args:
        dset (Dataset):  Model data.

    Returns:
        Dict: One temperature function for each station.
    """
    temperature_funcs = dict()
    for ivsname in dset.unique("ivsname"):
        try:
            site_id = dset.meta["station"][ivsname]["site_id"]
        except KeyError:
            site_id = ""

        idx_1 = dset.filter(ivsname_1=ivsname)
        time_1 = dset.time.utc[idx_1].mjd
        temp_1 = dset.temperature_1[idx_1]

        idx_2 = dset.filter(ivsname_2=ivsname)
        time_2 = dset.time.utc[idx_2].mjd
        temp_2 = dset.temperature_2[idx_2]

        time = np.append(time_1, time_2)
        temp = np.append(temp_1, temp_2)

        if any(np.isnan(temp)):
            vmf1_sta = apriori.get("vmf1_station", time=dset.time)
            try:
                temperature_funcs[ivsname] = vmf1_sta[ivsname]["temp"]
                log.warn(f"Missing temperature data for {ivsname} ({site_id}). Getting temperature from VMF1 station files")
                continue
            except KeyError:
                log.warn(f"Missing temperature data for {ivsname} ({site_id}). No backup data source available")

        # Amplitude, phase, offset, trend
        guess = [3 * np.std(temp) / (2 ** 0.5), 0, np.mean(temp), 0]
        optimize_func = lambda x: (x[0] * np.sin(2 * np.pi * time + x[1]) + x[2] + x[3] * time - temp)

        if len(temp) >= 4:
            coeff = optimize.leastsq(optimize_func, guess)[0]
        else:
            log.warn(f"Few datapoints, applying constant temperature for {ivsname} ({site_id}).")
            coeff = [0, 0, np.mean(temp), 0]
        temperature_funcs[ivsname] = _get_temperature_func(*coeff)
    return temperature_funcs


def _get_temperature_func(amplitude, phase, offset, trend):
    def temperature_func(time):
        return amplitude * np.sin(2 * np.pi * time + phase) + offset + trend * time

    return temperature_func
