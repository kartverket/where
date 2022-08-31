"""Calculate the delay caused by the troposphere for radio techniques

Description:
------------

Different models are available to determine the troposphere delay for radio techniques. The main models are GPT/GMF,
GPT2, GPT2w and VMF1 gridded.

References:
-----------
.. [1] IERS Conventions (2010) Software Library.
       http://maia.usno.navy.mil/conv2010/software.html



TODO:
    * Potential idea: Troposphere class. Remove meteorological data as a function. Create properties of pressure, 
      temperature, e, tm, lambd instead. Lazy evaluation based on config or input parameters. 
    * Test Dataset needed to assure better unit testing. So far all models are tested against output of the GIPSY
      program 'tropnominal', except GPT2w model.


"""
# External library imports
import numpy as np
import os

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where import apriori
from where.ext import iers_2010 as iers
from where.ext import gpt2w as ext_gpt2w
from where.lib import config
from where.lib import log

# Name of model
MODEL = __name__.split(".")[-1]

# Available troposphere models
MAPPING_FUNCTIONS = ["gmf", "gpt2", "gpt2w", "vmf1_gridded", "vmf1_station"]
METEOROLOGICAL_MODELS = ["vmf1_gridded", "vmf1_station", "gpt", "gpt2", "gpt2w", "site_pressure", "default"]
ZENITH_WET_DELAY_MODELS = ["none", "askne", "davis", "saastamoinen", "vmf1_gridded", "vmf1_station"]
ZENITH_HYDROSTATIC_DELAY_MODELS = ["saastamoinen", "vmf1_gridded", "vmf1_station"]
GRADIENT_MODELS = ["none", "apg"]

# Cache for GPT2 model
_GPT2 = dict()

# Cache for GPT2W model
_GPT2W = dict()


@plugins.register
def troposphere_for_all_stations(dset):
    """Calculate tropospheric delay for all stations

    Depending on the used technique a Dataset can include one or more stations. In the case of VLBI a baseline between
    two stations are analysed. Here the relative tropospheric delay is needed between two stations and not the absolute
    one.  Therefore if a Dataset has two stations then the difference between the tropospheric delay is calculated by
    using ``multiplier``, which is '-1' for station 1 and '1' for station 2. For one or more than two stations in a
    Dataset the ``multiplier`` is equal to '1'.

    Args:
        dset (Dataset):     Model data.

    Returns:
        numpy.ndarray: Corrections in meters for each observation.
    """
    dset_out = np.zeros(dset.num_obs)
    for multiplier in dset.for_each_suffix("station"):
        dset_out += multiplier * troposphere(dset)

    return dset_out


def troposphere(dset):
    r"""Calculate tropospheric delay for a station

    The tropospheric delay :math:`\Delta T` is calculated using equation (9.12) in the IERS Conventions
    :cite:`iers2010`,

    .. math:: \Delta T = \Delta T_h^z \cdot m_h(e) + \Delta T_w \cdot m_w + \Delta T_{asym} ,

    with :math:`m_h` and :math:`m_w` denoting the hydrostatic and wet mapping functions, :math:`\Delta T_h^z` zenith
    hydrostatic delay, :math:`\Delta T_w^z` zenith wet delay and :math:`\Delta T_{asym}` asymmetric delay.

    Different models for the mapping functions, the zenith hydrostatic and wet delay, and the asymmetric delay
    determination are available and are defined in the configuration file.

    The different terms of the tropospheric delay are stored in the Dataset in a table called ``troposphere`` in fields
    prefixed by ``troposphere_``.

    Args:
        dset (Dataset):     Model data.

    Returns:
        numpy.ndarray:  Total tropospheric delay for each observation and one station in [m].
    """
    
    log.debug(f"Computing troposphere for station{dset.default_field_suffix}")

    latitude, longitude, height = dset.site_pos.pos.llh.T
    time = dset.time
    zenith_distance = dset.site_pos.zenith_distance
    azimuth = dset.site_pos.azimuth
    elevation = dset.site_pos.elevation
    stations = dset.station

    obs_pressure = meteo_from_dset(dset, "pressure")
    obs_temp = meteo_from_dset(dset, "temperature")
    # TODO: use different names?
    obs_e = meteo_from_dset(dset, "e")
    obs_tm = meteo_from_dset(dset, "tm")
    obs_lambd = meteo_from_dset(dset, "lambd")

    pressure, temperature, e, tm, lambd = meteorological_data(stations, latitude, longitude, height, 
                                                              time, obs_pressure, obs_temp, 
                                                              obs_e, obs_tm, obs_lambd)
    #import IPython; IPython.embed()
    mh, mw = mapping_function(latitude, longitude, height, time, zenith_distance)
    mg, gn, ge = gradient_model(latitude, longitude, azimuth, elevation)
    zhd = zenith_hydrostatic_delay(stations, latitude, longitude, height, time, pressure)
    zwd = zenith_wet_delay(stations, latitude, longitude, height, time, temperature, e, tm, lambd)

    # Line of sight delay
    dT = mh * zhd + mw * zwd + mg * (gn * np.cos(azimuth) + ge * np.sin(azimuth))

    terms_and_levels = dict(
        dT=("operational", "meter"),
        mh=("operational", "dimensionless"),
        zhd=("operational", "meter"),
        mw=("operational", "dimensionless"),
        zwd=("operational", "meter"),
        mg=("operational", "dimensionless"),
        gn=("operational", "meter"),
        ge=("operational", "meter"),
    )
    for term, (write_level, unit) in terms_and_levels.items():
        field = "troposphere_{}{}".format(term, (dset.default_field_suffix or ""))
        if field in dset.fields:
            dset[field][:] = locals()[term]
        else:
            dset.add_float(field, val=locals()[term], write_level=write_level, unit=unit)

    # +DEBUG
    # if True:
    #    refdate_mjd = date.dt_to_mjd(dset.ref_date)
    #    delta_mjd = dset.time_unit.seconds / 86400
    #    for obs, sat in enumerate(dset.satellite):
    #        mjd = dset.time.utc.mjd[obs]
    #        zd = dset.site_pos.zenith_distance[obs]
    #        sat_id = dset.satellite_list[sat]
    #        print('MURKS: '+sat_id,mjd,zd,dT[obs],zhd[obs],zwd[obs],mh[obs],mw[obs],trop_asym[obs],pressure[obs])
    # -DEBUG

    return dT


def meteorological_data(stations, latitude, longitude, height, time,
                        obs_pressure=None,
                        obs_temp=None,
                        obs_e=None, 
                        obs_tm=None,
                        obs_lambd=None):
    """Determine meteorological data (atmospheric pressure, ...) based on configuration file definition.

    Args:
        stations (numpy.ndarray):        Station name for each observation.
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation
        obs_pressure (numpy.ndarray):    Observed pressure for each observation in [hPa]
        obs_temperature (numpy.ndarray): Observed temperature for each observation in [C]
        obs_e (numpy.ndarray):           Observed water vapor pressure for each observation in [hPa]
        obs_tm (numpy.ndarray):          Observed mean water vapor temperature for each observation in [K]
        obs_lambd (numpy.ndarray):       Observed water vapor decrease factor for each observation

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =====================================================================================
     Element       Unit         Description
    ============  ===========  =====================================================================================
     pressure      hPa          Atmospheric pressure value
     temperature   Celsius      Temperature values. Marked with 'None', if unknown.
     e             hPa          Water vapor pressure. Marked with 'None', if unknown.
     tm            K            Mean temperature of the water vapor. Marked with 'None', if unknown.
     lambd                      Water vapor decrease factor for each observation. Marked with 'None', if unknown.
    ============  ===========  =====================================================================================
    """

    models = config.tech[MODEL].meteorological_data.list
    if not models:
        log.fatal(
            f"No meteorological data model defined. "
            f"Available options are {', '.join(METEOROLOGICAL_MODELS)}"
        )
    log.debug(f"Meteorological data model: {models}")

    num_obs = len(time)
    
    pressure = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    temperature = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    e = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    tm = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    lambd = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    
    for model in models:
        if model == "vmf1_gridded":
            model_pressure = vmf1_gridded_pressure(latitude, longitude, height, time)
            model_temp, model_e, model_tm, model_lambd = None, None, None, None
        elif model == "vmf1_station":
            model_pressure = vmf1_station_pressure(stations, time)
            model_temp, model_e, model_tm, model_lambd = None, None, None, None
        elif model == "gpt":
            model_pressure, model_temp, _ = gpt(latitude, longitude, height, time)
            model_e, model_tm, model_lambd = None, None, None
        elif model == "gpt2":
            model_pressure, model_temp, _, model_e, _, _, _ = gpt2_meteo(latitude, longitude, height, time)
            model_tm, model_lambd = None, None, 
        elif model == "gpt2w":
            model_pressure, model_temp, _, model_tm, model_e, _, _, model_lambd, _ = gpt2w(latitude, longitude, height, time)
        elif model == "site_pressure":
            model_pressure = obs_pressure
            model_temp, model_e, model_tm, model_lambd = None, None, None, None
        elif model == "site_temperature":
            model_temp = obs_temperature
            model_pressure, model_e, model_tm, model_lambd = None, None, None, None
        elif model == "site_e":
            model_e = obs_e
            model_pressure, model_temp, model_tm, model_lambd = None, None, None, None
        elif model == "site_tm":
            model_tm = obs_tm
            model_pressure, model_temp, model_e, model_lambd = None, None, None, None
        elif model == "site_lambd":
            model_lambd = obs_lambd
            model_pressure, model_temp, model_e, model_tm = None, None, None, None
        elif model == "default":
            model_pressure = np.full(num_obs, fill_value=1013.25) # 1013.25 HPa (1 atm)
            # TODO: create default values for other parameters
            model_temp, model_e, model_tm, model_lambd = None, None, None, None
        else:
            log.fatal(
                f"Unknown meteorological data model '{model}'. Available models are {', '.join(METEOROLOGICAL_MODELS)}."
            )
        
        if model_pressure is not None:
            pressure = _update_data(pressure, model_pressure, model, "pressure")
        if model_temp is not None:
            temperature = _update_data(temperature, model_temp, model, "temperature")
        if model_e is not None:
            e = _update_data(e, model_e, model, "water vapor pressure")
        if model_tm is not None:
            tm = _update_data(tm, model_tm, model, "mean temperature of the water vapor")
        if model_lambd is not None:
            lambd = _update_data(lambd, model_lambd, model, "water vapor decrease factor")

    return pressure, temperature, e, tm, lambd


def gradient_model(latitude, longitude, azimuth, elevation):
    """Calculates asymmetric delay based on gradient model given in configuration file

    Args:
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):   Geodetic longitude for each observation in [rad]
        azimuth (numpy.ndarray):     Azimuth angle for each observatioin in [rad]
        elevation (numpy.ndarray):   Elevation angle for each observation in [rad]

    Returns:
        numpy.ndarray:  Troposphere asymmetric delay in [m] for each observation
    """
    models = config.tech[MODEL].gradients.list
    if not models:
        log.fatal(
            f"No gradient model defined. "
            f"Available options are {', '.join(GRADIENT_MODELS)}"
        )
    log.debug(f"Troposphere gradient model: {models}")

    
    # Note: Be aware, that the apg.f function uses c = 0.0031 instead of c = 0.0032.
    mg = 1 / (np.sin(elevation) * np.tan(elevation) + 0.0032)

    num_obs = len(latitude)
    gn = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    ge = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    
    for model in models:
        if model == "none":
            model_gn = np.zeros(num_obs)
            model_ge = np.zeros(num_obs)
        elif model == "apg":
            model_gn, model_ge = apg_gradient_model(latitude, longitude, azimuth, elevation)
        else:
            log.fatal(f"Unknown troposphere gradient model {model}. Available models are {', '.join(GRADIENT_MODELS)}")
    
        gn = _update_data(gn, model_gn, model, "north gradient")
        ge = _update_data(ge, model_ge, model, "east gradient")

    log.debug(f"Troposphere gradients North and East (average): {np.mean(gn)} {np.mean(ge)} [m]")
    return mg, gn, ge


def mapping_function(latitude, longitude, height, time, zenith_distance):
    """Calculates hydrostatic and wet mapping functions based on configuration file definition

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation
        zenith_distance (numpy.ndarray): Zenith distance for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =====================================================================================
     Element       Unit         Description
    ============  ===========  =====================================================================================
     mh                         Hydrostatic mapping function values
     mw                         Wet mapping function values
    ============  ===========  =====================================================================================
    """
    models = config.tech[MODEL].mapping_function.list
    if not models:
        log.fatal(
            f"No troposphere mapping functions defined. "
            f"Available options are {', '.join(MAPPING_FUNCTIONS)}"
        )
    log.debug(f"Troposphere mapping function: {models}")

    num_obs = len(time)
    mh = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    mw = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)
    
    for model in models:
        if model == "gmf":
            model_mh, model_mw = gmf_mapping_function(latitude, longitude, height, time, zenith_distance)
        elif model == "gpt2":
            _, _, _, _, model_mh, model_mw, _ = gpt2_mapping_function(latitude, longitude, height, time, zenith_distance)
        elif model == "gpt2w":
            _, _, _, _, _, model_mh, model_mw, _, _ = gpt2w_mapping_function(latitude, longitude, height, time, zenith_distance)
        elif model == "vmf1_gridded":
            model_mh, model_mw = vmf1_gridded_mapping_function(latitude, longitude, height, time, zenith_distance)
        elif model == "vmf1_station":
            model_mh, model_mw = vmf1_station_mapping_function(latitude, longitude, height, time, zenith_distance)
    
        else:
            log.fatal(
                f"Unknown troposphere mapping function {model}. "
                f"Available mapping functions are {', '.join(MAPPING_FUNCTIONS)}"
            )
        mh = _update_data(mh, model_mh, model, "hydrostatic mapping function")
        mw = _update_data(mw, model_mw, model, "wet mapping function")
    
    log.debug(f"Troposphere mapping functions mh, mw (average): {np.mean(mh)}, {np.mean(mw)}")
    return mh, mw


def zenith_hydrostatic_delay(stations, latitude, longitude, height, time, pressure):
    """Calculates zenith hydrostatic delay based on configuration file definition

    Args:
        stations (numpy.ndarray):    Station name for each observation.
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        ongitude (numpy.ndarray):    Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):      Orthometric height for each observation in [m]
        time (Time):                 Epoch of each observation
        pressure (numpy.ndarray):    Atmospheric pressure for each observation in [hPa]

    Returns:
        numpy.ndarray:         Array with zenith hydrostatic delay for each observation in [m]
    """

    models = config.tech[MODEL].zenith_hydrostatic_delay.list
    if not models:
        log.fatal(
            f"No zenith hydrostatic delay model defined. "
            f"Available options are {', '.join(ZENITH_HYDROSTATIC_DELAY_MODELS)}"
        )
    log.debug(f"Troposphere zenith hydrostatic delay model: {models}")
    
    num_obs = len(time)
    zhd = np.nan if num_obs == 1 else  np.full(num_obs, fill_value=np.nan)

    for model in models:
        if model == "saastamoinen":
            model_zhd = saastamoinen_zenith_hydrostatic_delay(pressure, latitude, height)
        elif model == "vmf1_gridded":
            model_zhd = vmf1_gridded_zenith_hydrostatic_delay(latitude, longitude, height, time)
        elif model == "vmf1_station":
            model_zhd = vmf1_station_zenith_hydrostatic_delay(stations, time)
        else:
            log.fatal(
                f"Unknown zenith hydrostatic troposphere delay model {model}. "
                f"Available models are {', '.join(ZENITH_HYDROSTATIC_DELAY_MODELS)}"
            )
        zhd = _update_data(zhd, model_zhd, model, "zenith hydrostatic delay")
        
    log.debug(f"Troposphere zenith hydrostatic delay (average): {np.mean(zhd)} [m]")
    return zhd


def zenith_wet_delay(stations, latitude, longitude, height, time, temperature, e, tm, lambd):
    """Calculates zenith wet delay based on configuration file definition

    Args:
        stations (numpy.ndarray):    Station name for each observation.
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):   Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):      Orthometric height for each observation in [m]
        time (Time):                 Epoch of each observation
        temperature (numpy.ndarray):   Temperature for each observation in [Celsius]
        e (numpy.ndarray):             Water vapor pressure for each observation in [hPa]
        tm (numpy.ndarray):            Mean temperature of the water vapor for each observation in [K]
        lambd (numpy.ndarray):         Water vapor decrease factor for each observation


    Returns:
        numpy.ndarray:    Zenith wet delay values for each observation in [m].
    """
    models = config.tech[MODEL].zenith_wet_delay.list
    if not models:
        log.fatal(
            f"No zenith wet delay model defined. "
            f"Available options are {', '.join(ZENITH_WET_DELAY_MODELS)}"
        )
    log.debug(f"Troposphere zenith wet delay model: {models}")

    num_obs = len(time)
    zwd = np.nan if num_obs == 1 else np.full(num_obs, fill_value=np.nan)

    for model in models:
        if model == "none":
            model_zwd = np.zeros(num_obs)
        elif model == "askne":
            model_zwd = askne_zenith_wet_delay(e, tm, lambd)
        elif model == "davis":
            model_zwd = davis_zenith_wet_delay(latitude, height, temperature, e)
        elif model == "saastamoinen":
            model_zwd = saastamoinen_zenith_wet_delay(latitude, height, temperature, e)
        elif model == "vmf1_gridded":
            model_zwd = vmf1_gridded_zenith_wet_delay(latitude, longitude, height, time)
        elif model == "vmf1_station":
            model_zwd = vmf1_station_zenith_wet_delay(stations, time)
        else:
            log.fatal(
                f"Unknown zenith wet troposphere delay model {model}. "
                f"Available models are {', '.join(ZENITH_WET_DELAY_MODELS)}"
            )

        zwd = _update_data(zwd, model_zwd, model, "zenith wet delay")

    log.debug(f"Troposphere zenith wet delay (average): {np.mean(zwd)} [m]")
    return zwd


def apg_gradient_model(latitude, longitude, azimuth, elevation):
    """Calculates ECMWF gradient model based on Fortran routine 'apg.f'

    The tropospheric asymmetric delay and the horizontal delay gradients ``G_N`` and ``G_E`` are computed by using the
    ``apg.f`` function from the IERS software library [1]_. We use the asymmetric delay, which is based on Equation
    (9.12) described in IERS conventions 2010 in Section 9.2 :cite:`iers2010`.

    Args:
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):   Geodetic longitude for each observation in [rad]
        azimuth (numpy.ndarray):     Azimuth angle for each observatioin in [rad]
        elevation (numpy.ndarray):   Elevation angle for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =====================================================================================
     Element       Unit         Description
    ============  ===========  =====================================================================================
     gn            m            Horizontal delay gradient in the North direction
     ge            m            Horizontal delay gradient in the East direction
    ============  ===========  =====================================================================================
    """
    num_obs = len(latitude)
    gn = np.empty(num_obs)
    ge = np.empty(num_obs)

    for obs in range(num_obs):
        # Get horizontal gradients G_N and G_E
        _, gn[obs], ge[obs] = iers.apg(latitude[obs], longitude[obs], azimuth[obs], elevation[obs])

    return gn * Unit.mm2m, ge * Unit.mm2m


def askne_zenith_wet_delay(e, tm, lambd):
    """Calculates zenith wet delay based on Askne and Nordius :cite:`askne1987`

    The Fortran routine 'asknewet.f' of the GPT2w library is used, which is based on Equation (18) of Askne and Nordius
    :cite:`askne1987`.

    Args:
        e (numpy.ndarray):      Water vapor pressure for each observation in [hPa]
        tm (numpy.ndarray):     Mean temperature of the water vapor for each observation in [K]
        lambd (numpy.ndarray):  Water vapor decrease factor for each observation

    Returns:
        numpy.ndarray:    Zenith wet delay for each observation in [m]
    """
    if not (len(e) == len(tm) == len(lambd)):
        log.fatal("Length of w, longitude and ellipsoidal height array is not equal.")

    num_obs = len(e)
    zwd = np.empty(num_obs)

    # Loop over all observations
    for obs in range(0, num_obs):
        zwd[obs] = ext_gpt2w.asknewet(e[obs], tm[obs], lambd[obs])

    return zwd


def davis_zenith_wet_delay(latitude, height, temperature, e):
    r"""Calculates zenith wet delay based on Saastamoinen/Davis model

    The total tropospheric delay for a given zenith distance :math:`z` is determined after Equation (19a) in
    Saastamoinen :cite:`saastamoinen1972`:

    .. math::
       \Delta T = 0.002277 \cdot \sec z \cdot (p + (1255/T + 0.05) \cdot e) - 1.16 \cdot \tan^2 z

    The zenith tropospheric delay is determined with :math:`z = 0`:

    .. math::
       \Delta T^z = \Delta T_h^z + \Delta T_w^z = zhd + zwd

    with the zenith hydrostatic delay :math:`zhd = 0.002277 \cdot p` and zenith wet delay :math:`zwd = 0.002277 \cdot
    (1255/T + 0.05) \cdot e`.

    A Fortran routine written by Davis corrects also for the gravity effect due to height :math:`H` and latitude
    :math:`\phi` (see http://acc.igs.org/tropo/wetsaas.f) and uses the constant :math:`0.0022768` instead of
    :math:`0.002277`, which leads to

    .. math::
       zwd = 0.0022768 \cdot (1255/T + 0.05) \cdot e / (1 - 0.00266 \cos 2 \phi - 0.00000028 H) .

    The difference between the orthometric height and geodetic height is called geoid undulation and can reach up to
    100 m due to Boehm et al. :cite:`boehm2007`. The influence of height difference is not significant in this
    equation.  The geodetic (ellipsoidal) height is therefore often used instead of the orthometric height.

    Args:
        latitude (numpy.ndarray):     Geodetic latitude for each observation in [rad]
        height (numpy.ndarray):       Orthometric height for each observation in [m]
        temperature (numpy.ndarray):  Temperature for each observation in [Celsius]
        e (numpy.ndarray):            Water vapor pressure for each observation in [hPa]

    Returns:
        numpy.ndarray:    Zenith wet delay for each observation in [m]
    """
    gravity_corr = saastamoinen_gravity_correction(latitude, height)

    # Zenith wet delay based on Eq. (19a) in Saastamoinen :cite:`saastamoinen1972` with additional gravity
    # correction
    zwd = 0.002_276_8 * (1255 / Unit.celsius_to_kelvin(temperature) + 0.05) * e / gravity_corr

    return zwd


def gmf_mapping_function(latitude, longitude, height, time, zenith_distance):
    """Calculates GMF hydrostatic and wet mapping functions

    Use the 'gmf.f' Fortran routine from the IERS software library to calculate the Global Mapping Function (see
    Section 9.2 in :cite:`iers2010`), which are described in Boehm et al. :cite:`boehm2006b`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation
        zenith_distance (numpy.ndarray): Zenith distance for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
    ============  ===========  =======================================================
    """
    num_obs = len(time)
    mh = np.empty(num_obs)
    mw = np.empty(num_obs)
    mjd = time.utc.mjd
    for obs in range(num_obs):
        mh[obs], mw[obs] = iers.gmf(mjd[obs], latitude[obs], longitude[obs], height[obs], zenith_distance[obs])

    return mh, mw


def gpt(latitude, longitude, height, time):
    """Calculates Global Pressure and Temperature (GPT)

    This subroutine determines atmospheric pressure and temperature globally based on spherical harmonics up to degree
    and order 9 (see Boehm et al. :cite:`boehm2007`) for a given latitude, longitude and ellipsoidal height.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     pressure      hPa          Pressure value
     temperature   Celsius      Temperature values
     geoid_undu    m            Geoid undulation (based on 9x9 EGM model)
    ============  ===========  =======================================================
    """
    num_obs = len(time)
    pressure = np.empty(num_obs)
    temperature = np.empty(num_obs)
    geoid_undu = np.empty(num_obs)

    # Note: For GPT2 and GPT2w linear interpolation is done between daily solutions. Here
    #      it does not seem to be necessary, because the performance of GPT Fortran routine
    #      is not so worse than for GPT2 and GPT2w.
    mjd = time.utc.mjd
    for obs in range(num_obs):
        pressure[obs], temperature[obs], geoid_undu[obs] = iers.gpt(mjd[obs], latitude[obs], longitude[obs], height[obs])

    return pressure, temperature, geoid_undu


def gpt2_meteo(latitude, longitude, height, time):
    """Calculates meteorological data based on GPT2 model

    The GPT2 model is described in Lagler et al. :cite:`lagler2013`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     pressure      hPa          Pressure value
     temperature   Celsius      Temperature values
     dt            degree/km    Temperature lapse rate
     e             hPa          Water vapor pressure
     geoid_undu    m            Geoid undulation (based on 9x9 EGM model)
    ============  ===========  =======================================================
    """
    num_obs = len(time)
    press = np.empty(num_obs)
    temp = np.empty(num_obs)
    dt = np.empty(num_obs)
    e = np.empty(num_obs)
    ah = np.empty(num_obs)
    aw = np.empty(num_obs)
    undu = np.empty(num_obs)

    # Determine GPT2 values for each observation by interpolating between two unique
    # daily solutions
    mjd = time.utc.mjd

    for obs in range(num_obs):
        # Start 'gpt2.f' day-by-day in folder where 'gpt2_5.grd' is placed and carry out
        # linear interpolation
        press[obs], temp[obs], dt[obs], e[obs], ah[obs], aw[obs], undu[obs] = gpt2_wrapper(
            mjd[obs], [latitude[obs]], [longitude[obs]], [height[obs]]
        )

    return press, temp, dt, e, undu

def gpt2_mapping_function(latitude, longitude, height, time, zenith_distance):
    """Calculates mapping function based on coefficients from GPT2 model

    The GPT2 model is described in Lagler et al. :cite:`lagler2013`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation
        zenith_distance (numpy.ndarray): Zenith distance for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
    ============  ===========  =======================================================
    """
    num_obs = len(time)
    ah = np.empty(num_obs)
    aw = np.empty(num_obs)
    mh = np.empty(num_obs)
    mw = np.empty(num_obs)

    # Determine GPT2 values for each observation by interpolating between two unique
    # daily solutions
    mjd = time.utc.mjd

    for obs in range(num_obs):
        # Start 'gpt2.f' day-by-day in folder where 'gpt2_5.grd' is placed and carry out
        # linear interpolation
        _, _, _, _, ah[obs], aw[obs], _ = gpt2_wrapper(
            mjd[obs], [latitude[obs]], [longitude[obs]], [height[obs]]
        )

        # Determine mapping function values based on coefficients 'ah' and 'aw'
        mh[obs], mw[obs] = iers.vmf1_ht(ah[obs], aw[obs], mjd[obs], latitude[obs], height[obs], zenith_distance[obs])

    return mh, mw

def gpt2_wrapper(mjd, latitude, longitude, hell):
    """Calculates meteorological data and mapping function coefficients based on GPT2 model

    The functions calls the IERS library routine 'gpt2.f' (see Section 9.2 in :cite:`iers2010`). The Fortran routine
    ``gpt2.f`` reads the grid file ``gpt2_5.grd``, which should be available in the same folder, where the Fortran
    programs runs. Therefore we change the current directory to the IERS source directory, so that 'gpt2.f' can read
    the grid file.

    Due to performance reasons the GPT2 values are not determined for each observation. The call of the Fortran routine
    ``gpt2.f`` takes time, because the grid file ``gpt2_5.grd`` has to be read for each observation. Instead the GPT2
    values are calculated only once for each unique day (modified Julian date rounded to integer) and saved in the
    cache _GPT2. The final GPT2 values are computed by a linear interpolation of the daily determined GPT2 values. The
    difference between the use of routine ``gpt2.f`` for each observation and the linear interpolation between daily
    solution is on the submillimeter level and can therefore be neglected.

    Args:
        mjd (numpy.float64):  Modified Julian date.
        latitude (list):      Array with latitude for each station in [rad].
        longitude (list):     Array with longitude for each station in [rad].
        hell (list):          Array with height for each station in [m].

    Returns:
        numpy.ndarray:  Array with following entries:

    =======  ===========  =======================================================
     Index    Unit         Description
    =======  ===========  =======================================================
     [0]      hPa          Pressure value
     [1]      Celsius      Temperature values
     [2]      degree/km    Temperature lapse rate
     [3]      hPa          Water vapor pressure
     [4]                   Hydrostatic mapping function coefficient ah
     [5]                   Wet mapping function coefficient aw
     [6]      m            Geoid undulation (based on 9x9 EGM model)
    =======  ===========  =======================================================
    """
    nstat = len(latitude)  # Number of stations
    it = 0  # Use of time variations (annual and semiannual terms)

    if not (len(latitude) == len(longitude) == len(hell)):
        log.fatal("Length of latitude, longitude and ellipsoidal height array is not equal.")


    # Change directory so that gpt2.f can read the gpt2_5.grd-file in the IERS source directory
    current_dir = os.getcwd()
    os.chdir(config.files.path(iers.__name__))

    # Loop over all unique dates (rounded to integer value)
    for date in _rounded_dates(mjd):

        # Check if date is already included in cache
        if date not in _GPT2:
            _GPT2[date] = np.array(iers.gpt2(date, latitude, longitude, hell, nstat, it)).reshape(-1)
    os.chdir(current_dir)

    # Linear interpolation between two daily GPT2 solutions
    mjd_int, mjd_frac = divmod(mjd, 1)
    output = _GPT2[mjd_int] + mjd_frac * (_GPT2[mjd_int + 1] - _GPT2[mjd_int])

    return output


def gpt2w_meteo(latitude, longitude, height, time):
    """Calculates meteorological data based on GPT2w model

    The GPT2w model is described in Boehm et al. :cite:`boehm2015`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     pressure      hPa          Pressure value
     temperature   Celsius      Temperature values
     dt            degree/km    Temperature lapse rate
     tm            K            Mean temperature of the water vapor
     e             hPa          Water vapor pressure
     la                         Water vapor decrease factor
     geoid_undu    m            Geoid undulation (based on 9x9 EGM model)
    ============  ===========  =======================================================
    """
    num_obs = len(time)
    press = np.empty(num_obs)
    temp = np.empty(num_obs)
    dt = np.empty(num_obs)
    tm = np.empty(num_obs)
    e = np.empty(num_obs)
    ah = np.empty(num_obs)
    aw = np.empty(num_obs)
    la = np.empty(num_obs)
    undu = np.empty(num_obs)

    mjd = time.utc.mjd

    # Determine GPT2W values for each observation by interpolating between two unique
    # daily solutions
    for obs in range(num_obs):
        # Start 'gpt2.f' day-by-day in folder where 'gpt2_5.grd' is placed and carry out
        # linear interpolation
        (press[obs], temp[obs], dt[obs], tm[obs], e[obs], ah[obs], aw[obs], la[obs], undu[obs]) = gpt2w_wrapper(
            mjd[obs], [latitude[obs]], [longitude[obs]], [height[obs]]
        )

        # Determine mapping function values based on coefficients 'ah' and 'aw'
        #mh[obs], mw[obs] = iers.vmf1_ht(ah[obs], aw[obs], mjd[obs], latitude[obs], height[obs], zenith_distance[obs])

    return press, temp, dt, tm, e, la, undu

def gpt2w_mapping_function(latitude, longitude, height, time, zenith_distance):
    """Calculates meteorological data and mapping function coefficients based on GPT2w model

    The GPT2w model is described in Boehm et al. :cite:`boehm2015`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation
        zenith_distance (numpy.ndarray): Zenith distance for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     pressure      hPa          Pressure value
     temperature   Celsius      Temperature values
     dt            degree/km    Temperature lapse rate
     tm            K            Mean temperature of the water vapor
     e             hPa          Water vapor pressure
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
     la                         Water vapor decrease factor
     geoid_undu    m            Geoid undulation (based on 9x9 EGM model)
    ============  ===========  =======================================================
    """
    num_obs = len(time)
    ah = np.empty(num_obs)
    aw = np.empty(num_obs)
    mh = np.empty(num_obs)
    mw = np.empty(num_obs)

    mjd = time.utc.mjd

    # Determine GPT2W values for each observation by interpolating between two unique
    # daily solutions
    for obs in range(num_obs):
        # Start 'gpt2.f' day-by-day in folder where 'gpt2_5.grd' is placed and carry out
        # linear interpolation
        _, _, _, _, _, ah[obs], aw[obs], _, _ = gpt2w_wrapper(
            mjd[obs], [latitude[obs]], [longitude[obs]], [height[obs]]
        )

        # Determine mapping function values based on coefficients 'ah' and 'aw'
        mh[obs], mw[obs] = iers.vmf1_ht(ah[obs], aw[obs], mjd[obs], latitude[obs], height[obs], zenith_distance[obs])

    return mh, mw



def gpt2w_wrapper(mjd, latitude, longitude, hell):
    """Calculates meteorological data and mapping function coefficients based on GPT2w model

    The functions calls the GPT2w library routine ``gpt2w_1w.f`` (see
    http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/GPT2w). The Fortran routine ``gpt2w_1w.f`` reads the grid file
    ``gpt2_1wA.grd``, which should be available in the same folder, where the Fortran programs runs. Therefore we
    change the current directory to the GPT2w source directory, so that ``gpt2w_1w.f`` can read the grid file.

    Due to performance reasons the GPT2w values are not determined for each observation. The call of the Fortran
    routine ``gpt2w_1w.f`` takes time, because the grid file ``gpt2_1wA.grd`` has to be read for each
    observation. Instead the GPT2w values are calculated only once for each unique day (modified Julian date rounded to
    integer) and saved in the cache _GPT2W. The final GPT2W values are computed by a linear interpolation of the daily
    determined GPT2W values.  The difference between the use of routine ``gpt2w_1w.f`` for each observation and the
    linear interpolation between daily solution is on the submillimeter level and can therefore be neglected.

    Args:
        mjd (numpy.float64):  Modified Julian date.
        latitude (list):      Array with latitude for each station in [rad].
        longitude (list):     Array with longitude for each station in [rad].
        hell (list):          Array with height for each station in [m].

    Returns:
        numpy.ndarray:  Array with following entries:

    =======  ===========  =======================================================
     Index    Unit         Description
    =======  ===========  =======================================================
     [0]      hPa          Pressure value
     [1]      Celsius      Temperature values
     [2]      degree/km    Temperature lapse rate
     [3]      K            Mean temperature of the water vapor
     [4]      hPa          Water vapor pressure
     [5]                   Hydrostatic mapping function coefficient ah
     [6]                   Wet mapping function coefficient aw
     [7]                   Water vapor decrease factor
     [8]      m            Geoid undulation (based on 9x9 EGM model)
    =======  ===========  =======================================================
    """
    nstat = len(latitude)  # Number of stations
    it = 0  # Use of time variations (annual and semiannual terms)

    if not (len(latitude) == len(longitude) == len(hell)):
        log.fatal("Length of latitude, longitute and ellipsoidal height array is not equal.")

    # Change directory so that gpt2w.f can read the gpt2_5.grd-file in the GPT2w source directory
    current_dir = os.getcwd()
    os.chdir(config.files.path(ext_gpt2w.__name__))

    # Loop over all unique dates (rounded to integer value)
    for date in _rounded_dates(mjd):

        # Check if date is already included in cache
        if date not in _GPT2W:
            _GPT2W[date] = np.array(ext_gpt2w.gpt2_1w(date, latitude, longitude, hell, nstat, it)).reshape(-1)
    os.chdir(current_dir)

    # Linear interpolation between two daily GPT2W solutions
    mjd_int, mjd_frac = divmod(mjd, 1)
    output = _GPT2W[mjd_int] + mjd_frac * (_GPT2W[mjd_int + 1] - _GPT2W[mjd_int])

    return output


def pressure_zhd(zhd, latitude, height):
    """Calculates atmospheric pressure based on zenith hydrostatic troposphere delay

    Args:
        zhd (numpy.ndarray):         Zenith hydrostatic delay for each observation in [m]
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        height (numpy.ndarray):      Orthometric height for each observation in [m]

    Returns:
        numpy.ndarray:    Atmospheric pressure for each observation in [hPa]
    """
    return saastamoinen_pressure(zhd, latitude, height)


def pressure_height_correction(pressure_ref, height_ref, height):
    r"""Atmospheric pressure is corrected from reference height to given height

    The atmospheric pressure height correction is based on a standard atmosphere model of Berg. See equation (1) in
    Boehm et al. :cite:`boehm2007` or equation (4) in Kouba :cite:`kouba2007`.

    The method uses equation (4) in Kouba :cite:`kouba2007` to determine the pressure for the reference height
    :math:`p_r` and for the station height :math:`p_s`:

    .. math::

       p_r &= 1013.25 (1 - 0.0000226 h_r)^{5.225} ,

       p_s &= 1013.25 (1 - 0.0000226 h_s)^{5.225} .

    The pressure at the station height is determined by setting these two equations
    in following relation:

    .. math:: p_s = (1 - 0.0000226h_r)^{5.225} / (1 - 0.0000226h_s)^{5.225} * p_r .

    TODO:
        The implementation of the pressure height correction is not completely described in Kouba :cite:`kouba2007`.  I
        found different solutions from Trimble, Per Helge and Bernese, whereby Bernese uses a similiar approach for the
        computation (see BERN52/LIB/FOR/VMF1ELL.f90). The use of the Trimble solution leads to differences of several
        millimeters for the zenith hydrostatic troposphere delay compared to our solution. Our solution compared to Per
        Helge's solution shows differences on the millimeter level.  What is the correct solution?

        Trimble solution::

            pressure = pressure_ref + 1013.25 * ((1 - 0.0000226 * height)**5.225 - (1 - 0.0000226 * height_ref)**5.225)

        Per Helge's solution::

            pressure = pressure_ref * (1 - 0.0000226 * (height-height_ref))**5.225

    Args:
        pressure_ref (numpy.ndarray):    Atmospheric pressure given for reference height for each observation in [hPa]
        height_ref (numpy.ndarray):      Reference height for each observation in [m]
        height (numpy.ndarray):          Height for each observation in [m]

    Returns:
        numpy.ndarray:      Pressure for given height for each observation in [hPa]
    """
    return ((1 - 0.000_022_6 * height) ** 5.225) / ((1 - 0.000_022_6 * height_ref) ** 5.225) * pressure_ref


def saastamoinen_gravity_correction(latitude, height):
    """Calculates gravity correction due to height and latitude

    The gravity correction due to height and latitude changes is shown in Davis et al. :cite:`davis1985` (see Equation
    (A13)) and is based on Saastamoinen :cite:`saastamoinen1972`.

    Args:
        latitude (numpy.ndarray):      Geodetic latitude for each observation in [rad]
        height (numpy.ndarray):        Orthometric height for each observation in [m]

    Returns:
        numpy.ndarray: Gravity correction for each observation in [m]
    """
    # Gravity correction based on Eq. (A13) in Davis et al. :cite:`davis1985`
    return 1 - 0.00266 * np.cos(2 * latitude) - 0.000_000_28 * height


def saastamoinen_pressure(zhd, latitude, height):
    """Calculates atmospheric pressure based on Saastamoinen model

    The difference between the orthometric height and geodetic height is called geoid undulation and can reach up to
    100 m due to Boehm et al. :cite:`boehm2007`. The influence of height difference is not significant by using
    equation (9.11). The geodetic (ellipsoidal) height is therefore often used instead of the orthometric height.

    Args:
        zhd (numpy.ndarray):         Zenith hydrostatic delay for each observation in [m]
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        height (numpy.ndarray):      Orthometric height for each observation in [m]

    Returns:
        numpy.ndarray:    Atmospheric pressure for each observation in [hPa]

    Example:
        >>> saastamoinen_pressure(2.2762945349647778, 0.83776, 200)
        1000.0
    """
    return zhd * saastamoinen_gravity_correction(latitude, height) / 0.002_276_8


def saastamoinen_zenith_hydrostatic_delay(pressure, latitude, height):
    r"""Calculates zenith hydrostatic delay with Saastamoinen model

    The zenith hydrostatic delay :math:`zhd` is calculated using equation (9.11) from Saastamoinen in the IERS
    Conventions :cite:`iers2010`:

    .. math::
         \Delta T_h^z = zhd = (0.0022768 \pm 0.0000005) P_0 / f_s(\phi, H)

    with the atmospheric pressure :math:`P_0` and the correction :math:`f_s(\phi, H)`
    (see equation (9.4) in :cite:`iers2010`):

    .. math::
       f_s(\phi, H) = 1 - 0.00266 \cos 2 \phi - 0.00000028 H ,

    where :math:`\phi` is the geodetic latitude of the station and :math:`H` is the geodetic height of the station.

    The difference between the orthometric height and geodetic height is called geoid undulation and can reach up to
    100 m due to Boehm et al. :cite:`boehm2007`. The influence of height difference is not significant by using
    equation (9.11). The geodetic (ellipsoidal) height is therefore often used instead of the orthometric height.

    Args:
        pressure (numpy.ndarray):      Atmospheric pressure for each observation in [hPa]
        latitude (numpy.ndarray):      Geodetic latitude for each observation in [rad]
        height (numpy.ndarray):        Orthometric height for each observation in [m]

    Returns:
        numpy.ndarray:    Zenith hydrostatic delay for each observation in [m]

    Example:
        >>> saastamoinen_zenith_hydrostatic_delay(1000, 0.83776, 200)
        2.2762945349647778

    TODO:
        Example is not compatible to http://ggosatm.hg.tuwien.ac.at/DELAY/SOURCE/GPT2w/FORTRAN/saasthyd.f
        example (zhd = 2.2695 m).
    """
    return 0.002_276_8 * pressure / saastamoinen_gravity_correction(latitude, height)


def saastamoinen_zenith_wet_delay(latitude, height, temperature, e):
    r"""Calculates zenith wet delay based Saastamoinen model

    The total tropospheric delay for a given zenith distance :math:`z` is determined after Equation (19a) in
    Saastamoinen :cite:`saastamoinen1972`:

    .. math::
       \Delta T = 0.002277 \cdot \sec z \cdot (p + (1255/T + 0.05) \cdot e) - 1.16 \cdot \tan^2 z

    The zenith tropospheric delay is determined with :math:`z = 0`:

    .. math::
       \Delta T^z = \Delta T_h^z + \Delta T_w^z = zhd + zwd

    with the zenith hydrostatic delay :math:`zhd = 0.002277 \cdot p` and zenith wet delay
    :math:`zwd = 0.002277 \cdot (1255/T + 0.05) \cdot e`.

    Args:
        latitude (numpy.ndarray):     Geodetic latitude for each observation in [rad]
        height (numpy.ndarray):       Orthometric height for each observation in [m]
        temperature (numpy.ndarray):  Temperature for each observation in [Celsius]
        e (numpy.ndarray):            Water vapor pressure for each observation in [hPa]

    Returns:
        numpy.ndarray:     Zenith wet delay for each observation in [m]
    """
    # Zenith wet delay based on Eq. (19a) in Saastamoinen :cite:`saastamoinen1972`
    return 0.002_276_8 * (1255 / Unit.celsius_to_kelvin(temperature) + 0.05) * e


def meteo_from_dset(dset, fieldname):
    """Get atmospheric pressure from local site measurements

    Missing values are set to np.nan

    Args:
        dset (Dataset): A Dataset containing model data.

    Returns:
        numpy.ndarray: Atmospheric pressure for each observation in [hPa]
    """
    meteo_data = np.full(dset.num_obs, fill_value=np.nan)

    if fieldname + (dset.default_field_suffix or "") in dset.fields:
        meteo_data = dset[fieldname]

    return meteo_data


def vmf1_gridded_pressure(latitude, longitude, height, time):
    """Calulates VMF1 gridded atmospheric pressure

    Equation (9.11) in IERS conventions 2010 :cite:`iers2010` is converted to 'pressure' and the VMF1 gridded
    atmospheric pressure can be determined with help of the gridded zenith hydrostatic delays from VMF1. The gridded
    VMF1 pressure values needs to be rescaled to actual station height. The method is described with Eq. (4) in Kouba
    :cite:`kouba2007`.

    Args:
        latitude (numpy.ndarray):    Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):   Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):      Orthometric height for each observation in [m]
        time (Time):                 Epoch of each observation

    Returns:
        numpy.ndarray:  Atmospheric pressure for each observation in [hPa].
    """
    # Get gridded VMF1 data
    vmf1 = apriori.get("vmf1_grid", time=time)
    try:
        grid_zhd = vmf1["zh"](time, longitude, latitude)  # Interpolation in time and space in VMF1 grid
        grid_height = vmf1["ell"](longitude, latitude, grid=False)
        grid_pressure = pressure_zhd(grid_zhd, latitude, grid_height)

        # Rescale gridded pressure to station pressure
        pressure = pressure_height_correction(grid_pressure, grid_height, height)
    except KeyError:
        pressure = np.full(len(time), fill_value=np.nan)
    return pressure

def vmf1_station_pressure(stations, time):
    """Calculates atmospheric pressure from station dependent VMF1 files

    Values for unknown station names are set to np.nan.

    Args:
        stations (numpy.ndarray):    Station name for each observation.
        time (Time):                 Epoch of each observation

    Returns:
        numpy.ndarray:  Atmospheric pressure for each observation in [hPa].
    """
    
    vmf1 = apriori.get("vmf1_station", time=time)
    
    if time.size == 1:
        try:
            return float(vmf1[stations]["pressure"](time.mjd))
        except KeyError:
            return np.nan
    
    pressure = np.full(len(time), fill_value=np.nan)

    for sta in np.unique(stations):
        idx = stations == sta
        try:
            pressure[idx] = vmf1[sta]["pressure"](time.mjd[idx])
        except KeyError:
            # pressure[idx] is already set to nan
            pass
    return pressure

def vmf1_gridded_mapping_function(latitude, longitude, height, time, zenith_distance):
    """Calculates VMF1 hydrostatic and wet mapping functions based on gridded VMF1 files

    This routine determines the VMF1 (Vienna Mapping Functions 1) described in Boehm et al. :cite:`boehm2006a` and
    Kouba :cite:`kouba2007` and uses the 'vmf1_ht.f' from the IERS software library :cite:`iers2010`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        longitude (numpy.ndarray):       Geodetic longitude for each observation in [rad]
        height (numpy.ndarray):          Orthometric height for each observation in [m]
        time (Time):                     Epoch of each observation
        zenith_distance (numpy.ndarray): Zenith distance for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
    ============  ===========  =======================================================
    """
    # Get gridded VMF1 data
    vmf1 = apriori.get("vmf1_grid", time=time)

    num_obs = len(time)
    mh = np.full(num_obs, fill_value=np.nan)
    mw = np.full(num_obs, fill_value=np.nan)

    for obs in range(num_obs):
        try:
            mh[obs], mw[obs] = iers.vmf1_ht(
                vmf1["ah"](time[obs], longitude[obs], latitude[obs]),
                vmf1["aw"](time[obs], longitude[obs], latitude[obs]),
                time.utc.mjd_int[obs],
                latitude[obs],
                height[obs],
                zenith_distance[obs],
            )
        except KeyError:
            # mh[obs] already set to nan
            pass

    return mh, mw

def vmf1_station_mapping_function(latitude, stations, time, zenith_distance):
    """Calculates VMF1 hydrostatic and wet mapping functions based on station VMF1 files

    This routine determines the VMF1 (Vienna Mapping Functions 1) described in Boehm et al. :cite:`boehm2006a` and
    Kouba :cite:`kouba2007` and uses the 'vmf1_ht.f' from the IERS software library :cite:`iers2010`.

    Args:
        latitude (numpy.ndarray):        Geodetic latitude for each observation in [rad]
        stations (numpy.ndarray):        Station name for each observation
        time (Time):                     Epoch of each observation
        zenith_distance (numpy.ndarray): Zenith distance for each observation in [rad]

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
    ============  ===========  =======================================================
    """
    vmf1 = apriori.get("vmf1_station", time=time)

    num_obs = len(time)
    mh = np.full(num_obs, fill_value=np.nan)
    mw = np.full(num_obs, fill_value=np.nan)
    mjd = time.utc.mjd

    for obs in range(num_obs):
        try:
            mh[obs], mw[obs] = iers.vmf1(
                vmf1[stations[obs]]["ah"](mjd[obs]),
                vmf1[stations[obs]]["aw"](mjd[obs]),
                mjd[obs],
                latitude[obs],
                zenith_distance[obs],
            )
        except KeyError:
            # mh[obs] already set to nan
            pass

    return mh, mw


def vmf1_gridded_zenith_wet_delay(latitude, longitude, height, time):
    """Calculates zenith wet delay based on gridded zenith wet delays from VMF1

    Uses gridded zenith wet delays from VMF1, which are rescaled from the gridded height to actual station height by
    using Equation(5) described in Kouba :cite:`kouba2007`.

    Args:
        latitude (numpy.ndarray):    Latitude for each observation
        longitude (numpy.ndarray):    Longitude for each observation
        height (numpy.ndarray): Height for each observation
        time (Time):            Epoch of each observation

    Returns:
        numpy.ndarray:     Zenith wet delay for each observation in [m]
    """
    # Get gridded VMF1 data
    vmf1 = apriori.get("vmf1_grid", time=time)
    
    try:
        grid_zwd = vmf1["zw"](time, longitude, latitude)  # Interpolation in time and space in VMF1 grid
        grid_height = vmf1["ell"](longitude, latitude, grid=False)

        # Zenith Wet delay. Eq. (5) in Kouba :cite:`kouba2007`
        zwd = grid_zwd * np.exp(-(height - grid_height) / 2000)
    except KeyError:
        zwd = np.full(num_obs, fill_value=np.nan)

    return zwd

def vmf1_station_zenith_wet_delay(stations, time):
    """Calculates zenith wet delay based on station zenith wet delays from VMF1

    Args:
        stations (numpy.ndarray):   Station name for each observation
        time (Time):                Epoch of each observation

    Returns:
        numpy.ndarray:     Zenith wet delay for each observation in [m]
    """
    vmf1 = apriori.get("vmf1_station", time=time)
    zwd = np.full(len(time), fill_value=np.nan)

    for sta in np.unique(stations):
        idx = stations == sta
        try:
            zwd[idx] = vmf1[sta]["zw"](time.mjd[idx])
        except KeyError:
            # zwd[obs] already set to nan
            pass
    return zwd

def vmf1_gridded_zenith_hydrostatic_delay(latitude, longitude, height, time):
    """Calculates zenith hydrostatic delay based on gridded zenith wet delays from VMF1

    Uses gridded zenith wet delays from VMF1, which are rescaled from the gridded height to actual station height by
    using Equation (3) and (4) described in Kouba :cite:`kouba2007`.

    Args:
        latitude (numpy.ndarray):    Latitude for each observation
        longitude (numpy.ndarray):    Longitude for each observation
        height (numpy.ndarray): Height for each observation
        time (Time):            Epoch of each observation

    Returns:
        numpy.ndarray:     Zenith hydrostatic delay for each observation in [m]
    """
    # Gridded pressure rescaled to station height
    pressure = vmf1_gridded_pressure(latitude, longitude, height, time)
    
    # Zenith hydrostatic delay. Eq. (4) in Kouba :cite:`kouba2007`
    zhd = saastamoinen_zenith_hydrostatic_delay(pressure, latitude, height)

    return zhd

def vmf1_station_zenith_hydrostatic_delay(stations, time):
    """Calculates zenith hydrostatic delay based on station zenith hydrostatic delays from VMF1

    Args:
        stations (numpy.ndarray):   Station name for each observation.
        time (Time):                Epoch of each observation

    Returns:
        numpy.ndarray:     Zenith hydrostatic delay for each observation in [m]
    """
    vmf1 = apriori.get("vmf1_station", time=time)
    zhd = np.full(len(time), fill_value=np.nan)

    for sta in np.unique(stations):
        idx = stations == sta
        try:
            zhd[idx] = vmf1[sta]["zh"](time.mjd[idx])
        except KeyError:
            # zhd[idx] already set to nan
            pass

    return zhd

def _update_data(data, model_data, model, model_text):
    """Updates the input argument data with information from model_data if entries are missing in data.
    
    Args:
        data (numpy.ndarray or float):        Data accumulated so far
        model_data(numpy.ndarray or float):   New data
        model (str):                          Name of model
        model_text (str):                     Text description of model
    
    Returns:
        np.ndarray or float:                  Data updated with new data
    """
    idx_missing = np.isnan(data)
    num_missing = np.sum(idx_missing)
    if num_missing > 0 and isinstance(model_data, np.ndarray):
        data[idx_missing] = model_data[idx_missing]
        num_missing_model = np.sum(np.isnan(model_data))
        log.debug(f"Used {model} for {model_text} for {num_missing-num_missing_model} observations")
    elif num_missing > 0 and isinstance(model_data, float):
        data = model_data
        log.debug(f"Used {model} for {model_text}")

    return data

def _rounded_dates(datetimes):
    """Calculate days including the given datetimes

    Args:
        datetimes (numpy.float64):   np.array with datetimes as floats

    Returns:
        numpy.ndarray:  Unique days representing the datetimes.

    Example:
        >>> _rounded_dates(np.array([54001.2, 54001.4, 54001.8, 54002, 54005.5, 54006.]))
        array([ 54001.,  54002.,  54003.,  54005.,  54006.,  54007.])
    """
    return np.unique(np.hstack((np.floor(datetimes), np.floor(datetimes) + 1)))
