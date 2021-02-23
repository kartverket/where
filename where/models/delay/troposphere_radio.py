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
    * Maybe it could be better to implement a "troposphere" class. For example with an object the information about
      used models can be retrieved, which is useful for example in routine :func:`zenith_wet_delay` (see TODO).
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
METEOROLOGICAL_MODELS = ["vmf1_gridded", "vmf1_station", "gpt", "gpt2", "gpt2w", "site_pressure"]
ZENITH_WET_DELAY_MODELS = ["none", "askne", "davis", "saastamoinen", "vmf1_gridded", "vmf1_station"]
ZENITH_HYDROSTATIC_MODELS = ["saastamoinen", "vmf1_gridded", "vmf1_station"]
GRADIENT_MODELS = ["none", "apg"]

# Default relation between used mapping function model and zenith wet delay model
MAPPING_ZENITH_WET_RELATION = dict(gmf="none", gpt2="saastamoinen", gpt2w="askne", vmf1_gridded="vmf1_gridded")

# Default relation between used mapping function model and meteorological data
MAPPING_METEO_RELATION = dict(gmf="gpt", gpt2="gpt2", gpt2w="gpt2w", vmf1_gridded="vmf1_gridded")

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
    latitude, _, height = dset.site_pos.pos.llh.T
    pressure, temperature, e, tm, lambd = meteorological_data(dset)
    mh, mw = mapping_function(dset)
    mg, gn, ge = gradient_model(dset)
    zhd = zenith_hydrostatic_delay(dset, pressure, latitude, height)
    zwd = zenith_wet_delay(dset, temperature, e, tm, lambd)

    # Line of sight delay
    dT = mh * zhd + mw * zwd + mg * (gn * np.cos(dset.site_pos.azimuth) + ge * np.sin(dset.site_pos.azimuth))

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
            # dset.add_float(field, table="troposphere", val=locals()[term], write_level=write_level)
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


def meteorological_data(dset):
    """Determine meteorological data (atmospheric pressure, ...) based on configuration file definition.

    Args:
        dset (Dataset): Model data.

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
    pressure = None
    temperature = None
    e = None
    tm = None
    lambd = None

    model = config.tech.get("meteorological_data", section=MODEL, default="").str
    mapping_function = config.tech[MODEL].mapping_function.str

    # Use default meteorological models, if no model is defined in configuration file
    if not model:
        try:
            model = MAPPING_METEO_RELATION[mapping_function]
        except KeyError:
            log.fatal(
                f"Unknown mapping function {mapping_function}. "
                f"Available mapping functions are {', '.join(MAPPING_FUNCTIONS)}"
            )

    log.debug(f"Meteorological data model: {model}")

    if model == "vmf1_gridded":
        pressure = vmf1_gridded_pressure(dset)
    if model == "vmf1_station":
        pressure = vmf1_station_pressure(dset)
    elif model == "gpt":
        pressure, temperature, _ = gpt(dset)

    elif model == "gpt2":
        pressure, temperature, _, e, _, _, _ = gpt2(dset)

    elif model == "gpt2w":
        pressure, temperature, _, tm, e, _, _, lambd, _ = gpt2w(dset)

    elif model == "site_pressure":
        pressure = site_pressure(dset)

    else:
        log.fatal(
            f"Unknown meteorological data model {model}. Available models are {', '.join(METEOROLOGICAL_MODELS)}"
        )

    return pressure, temperature, e, tm, lambd


def gradient_model(dset):
    """Calculates asymmetric delay based on gradient model given in configuration file

    Args:
        dset (Dataset):       Model data.

    Returns:
        numpy.ndarray:  Troposphere asymmetric delay in [m] for each observation
    """
    model = config.tech.get("gradients", section=MODEL, default="apg").str
    log.debug(f"Troposphere gradient model: {model}")

    # Note: Be aware, that the apg.f function uses c = 0.0031 instead of c = 0.0032.
    mg = 1 / (np.sin(dset.site_pos.elevation) * np.tan(dset.site_pos.elevation) + 0.0032)

    if model == "none":
        gn = ge = np.zeros(dset.num_obs)
    elif model == "apg":
        gn, ge = apg_gradient_model(dset)
    else:
        log.fatal(f"Unknown troposphere gradient model {model}. Available models are {', '.join(GRADIENT_MODELS)}")

    log.debug(f"Troposphere gradients North and East (average): {np.mean(gn)} {np.mean(ge)} [m]")

    return mg, gn, ge


def mapping_function(dset):
    """Calculates hydrostatic and wet mapping functions based on configuration file definition

    Args:
        dset (Dataset):    A Dataset containing model data.

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =====================================================================================
     Element       Unit         Description
    ============  ===========  =====================================================================================
     mh                         Hydrostatic mapping function values
     mw                         Wet mapping function values
    ============  ===========  =====================================================================================
    """
    model = config.tech[MODEL].mapping_function.str
    log.debug(f"Troposphere mapping function: {model}")

    if model == "gmf":
        mh, mw = gmf_mapping_function(dset)

    elif model == "gpt2":
        _, _, _, _, mh, mw, _ = gpt2(dset)

    elif model == "gpt2w":
        _, _, _, _, _, mh, mw, _, _ = gpt2w(dset)

    elif model == "vmf1_gridded":
        mh, mw = vmf1_gridded_mapping_function(dset)
    elif model == "vmf1_station":
        mh, mw = vmf1_station_mapping_function(dset)

    else:
        log.fatal(
            f"Unknown troposphere mapping function {model}. "
            f"Available mapping functions are {', '.join(MAPPING_FUNCTIONS)}"
        )

    return mh, mw


def zenith_hydrostatic_delay(dset, pressure, latitude, height):
    """Calculates zenith hydrostatic delay based on configuration file definition

    Args:
        pressure:    Array with atmospheric pressure for each observation in [hPa]
        latitude:    Array with geodetic latitude for each observation in [rad]
        height:      Array with orthometric height for each observation in [m]

    Returns:
        numpy.ndarray:         Array with zenith hydrostatic delay for each observation in [m]
    """
    model = config.tech.get("zenith_hydrostatic_delay", section=MODEL, default="saastamoinen").str
    log.debug(f"Troposphere zenith hydrostatic delay model: {model}")

    if model == "saastamoinen":
        zhd = saastamoinen_zenith_hydrostatic_delay(pressure, latitude, height)
    elif model == "vmf1_gridded":
        zhd = vmf1_gridded_zenith_hydrostatic_delay(dset)
    elif model == "vmf1_station":
        zhd = vmf1_station_zenith_hydrostatic_delay(dset)
    else:
        log.fatal(
            f"Unknown zenith hydrostatic troposphere delay model {model}. "
            f"Available models are {', '.join(ZENITH_HYDROSTATIC_MODELS)}"
        )

    log.debug(f"Troposphere zenith hydrostatic delay (average): {np.mean(zhd)} [m]")

    return zhd


def zenith_wet_delay(dset, temperature, e, tm, lambd):
    """Calculates zenith wet delay based on configuration file definition

    Args:
        dset (Dataset):                A Dataset containing model data
        temperature (numpy.ndarray):   Temperature for each observation in [Celsius]
        e (numpy.ndarray):             Water vapor pressure for each observation in [hPa]
        tm (numpy.ndarray):            Mean temperature of the water vapor for each observation in [K]
        lambd (numpy.ndarray):         Water vapor decrease factor for each observation


    Returns:
        numpy.ndarray:    Zenith wet delay values for each observation in [m].
    """
    model = config.tech.get("zenith_wet_delay", section=MODEL, default="").str
    mapping_function = config.tech[MODEL].mapping_function.str

    # Use default zenith wet delay models, if no model is defined in configuration file
    if not model:
        try:
            model = MAPPING_ZENITH_WET_RELATION[mapping_function]
        except KeyError:
            log.fatal(
                f"Unknown mapping function {mapping_function}. "
                f"Available mapping functions are {', '.join(MAPPING_FUNCTIONS)}"
            )
    log.debug(f"Troposphere zenith wet delay model: {model}")

    if model == "none":
        zwd = np.zeros(dset.num_obs)
    elif model == "askne":
        zwd = askne_zenith_wet_delay(e, tm, lambd)
        # TODO: log.fatal(f"Meteorological model {met_model!r} does not provide input parameters for using Askne and "
        #                 "Nordius zenith wet delay model. Use 'gpt2w' model")

    elif model == "davis":
        latitude, _, height = dset.site_pos.pos.llh.T
        zwd = davis_zenith_wet_delay(latitude, height, temperature, e)

    elif model == "saastamoinen":
        latitude, _, height = dset.site_pos.pos.llh.T
        zwd = saastamoinen_zenith_wet_delay(latitude, height, temperature, e)
        # TODO: log.fatal(f"Meteorological model {met_model} does not provide input parameters for using Saastamoinen "
        #                 "zenith wet delay model. Use 'gpt2' or 'gpt2w' model")

    elif model == "vmf1_gridded":
        zwd = vmf1_gridded_zenith_wet_delay(dset)
    elif model == "vmf1_station":
        zwd = vmf1_station_zenith_wet_delay(dset)
    else:
        log.fatal(
            f"Unknown zenith wet troposphere delay model {model}. "
            f"Available models are {', '.join(ZENITH_WET_DELAY_MODELS)}"
        )

    log.debug(f"Troposphere zenith wet delay (average): {np.mean(zwd)} [m]")

    return zwd


def apg_gradient_model(dset):
    """Calculates ECMWF gradient model based on Fortran routine 'apg.f'

    The tropospheric asymmetric delay and the horizontal delay gradients ``G_N`` and ``G_E`` are computed by using the
    ``apg.f`` function from the IERS software library [1]_. We use the asymmetric delay, which is based on Equation
    (9.12) described in IERS conventions 2010 in Section 9.2 :cite:`iers2010`.

    Args:
        dset (Dataset):       A Dataset containing model data.

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =====================================================================================
     Element       Unit         Description
    ============  ===========  =====================================================================================
     gn            m            Horizontal delay gradient in the North direction
     ge            m            Horizontal delay gradient in the East direction
    ============  ===========  =====================================================================================
    """
    gn = ge = np.empty(dset.num_obs)
    lat, lon, _ = dset.site_pos.pos.llh.T
    az = dset.site_pos.azimuth
    el = dset.site_pos.elevation

    for obs in range(dset.num_obs):
        # Get horizontal gradients G_N and G_E
        _, gn[obs], ge[obs] = iers.apg(lat[obs], lon[obs], az[obs], el[obs])

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


def gmf_mapping_function(dset):
    """Calculates GMF hydrostatic and wet mapping functions

    Use the 'gmf.f' Fortran routine from the IERS software library to calculate the Global Mapping Function (see
    Section 9.2 in :cite:`iers2010`), which are described in Boehm et al. :cite:`boehm2006b`.

    Args:
        dset (Dataset):    Model data.

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
    ============  ===========  =======================================================
    """
    mh = np.empty(dset.num_obs)
    mw = np.empty(dset.num_obs)
    lat, lon, height = dset.site_pos.pos.llh.T
    zd = dset.site_pos.zenith_distance
    mjd = dset.time.utc.mjd
    for obs in range(dset.num_obs):
        mh[obs], mw[obs] = iers.gmf(mjd[obs], lat[obs], lon[obs], height[obs], zd[obs])

    return mh, mw


def gpt(dset):
    """Calculates Global Pressure and Temperature (GPT)

    This subroutine determines atmospheric pressure and temperature globally based on spherical harmonics up to degree
    and order 9 (see Boehm et al. :cite:`boehm2007`) for a given latitude, longitude and ellipsoidal height.

    Args:
        dset (Dataset):    Model data.

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
    pressure = np.empty(dset.num_obs)
    temperature = np.empty(dset.num_obs)
    geoid_undu = np.empty(dset.num_obs)

    # Note: For GPT2 and GPT2w linear interpolation is done between daily solutions. Here
    #      it does not seem to be necessary, because the performance of GPT Fortran routine
    #      is not so worse than for GPT2 and GPT2w.
    mjd = dset.time.utc.mjd
    lat, lon, height = dset.site_pos.pos.llh.T
    for obs in range(dset.num_obs):
        pressure[obs], temperature[obs], geoid_undu[obs] = iers.gpt(mjd[obs], lat[obs], lon[obs], height[obs])

    return pressure, temperature, geoid_undu


def gpt2(dset):
    """Calculates meteorological data and mapping function coefficients based on GPT2 model

    The GPT2 model is described in Lagler et al. :cite:`lagler2013`.

    Args:
        dset (Dataset):    Model data.

    Returns:
        tuple of Numpy Arrays: Includes the following elements, each with entries for each observation

    ============  ===========  =======================================================
     Element       Unit         Description
    ============  ===========  =======================================================
     pressure      hPa          Pressure value
     temperature   Celsius      Temperature values
     dt            degree/km    Temperature lapse rate
     e             hPa          Water vapor pressure
     mh                         Hydrostatic mapping function coefficient ah
     mw                         Wet mapping function coefficient aw
     geoid_undu    m            Geoid undulation (based on 9x9 EGM model)
    ============  ===========  =======================================================
    """
    press = np.empty(dset.num_obs)
    temp = np.empty(dset.num_obs)
    dt = np.empty(dset.num_obs)
    e = np.empty(dset.num_obs)
    ah = np.empty(dset.num_obs)
    aw = np.empty(dset.num_obs)
    mh = np.empty(dset.num_obs)
    mw = np.empty(dset.num_obs)
    undu = np.empty(dset.num_obs)

    # Determine GPT2 values for each observation by interpolating between two unique
    # daily solutions
    mjd = dset.time.utc.mjd
    lat, lon, height = dset.site_pos.pos.llh.T
    zd = dset.site_pos.zenith_distance

    for obs in range(dset.num_obs):
        # Start 'gpt2.f' day-by-day in folder where 'gpt2_5.grd' is placed and carry out
        # linear interpolation
        press[obs], temp[obs], dt[obs], e[obs], ah[obs], aw[obs], undu[obs] = gpt2_wrapper(
            mjd[obs], [lat[obs]], [lon[obs]], [height[obs]]
        )

        # Determine mapping function values based on coefficients 'ah' and 'aw'
        mh[obs], mw[obs] = iers.vmf1_ht(ah[obs], aw[obs], mjd[obs], lat[obs], height[obs], zd[obs])

    return press, temp, dt, e, mh, mw, undu


def gpt2_wrapper(mjd, lat, lon, hell):
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
        lat (list):           Array with latitude for each station in [rad].
        lon (list):           Array with longitude for each station in [rad].
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
    nstat = len(lat)  # Number of stations
    it = 0  # Use of time variations (annual and semiannual terms)

    if not (len(lat) == len(lon) == len(hell)):
        log.fatal("Length of latitude, longitude and ellipsoidal height array is not equal.")

    # TODO: Case not handled if several stations are included in Dataset. Is that the case for VLBI?

    # Change directory so that gpt2.f can read the gpt2_5.grd-file in the IERS source directory
    current_dir = os.getcwd()
    os.chdir(config.files.path(iers.__name__))

    # Loop over all unique dates (rounded to integer value)
    for date in _rounded_dates(mjd):

        # Check if date is already included in cache
        if date not in _GPT2:
            _GPT2[date] = np.array(iers.gpt2(date, lat, lon, hell, nstat, it)).reshape(-1)
    os.chdir(current_dir)

    # Linear interpolation between two daily GPT2 solutions
    mjd_int, mjd_frac = divmod(mjd, 1)
    output = _GPT2[mjd_int] + mjd_frac * (_GPT2[mjd_int + 1] - _GPT2[mjd_int])

    return output


def gpt2w(dset):
    """Calculates meteorological data and mapping function coefficients based on GPT2w model

    The GPT2w model is described in Boehm et al. :cite:`boehm2015`.

    Args:
        dset (Dataset):    Model data.

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
    press = np.empty(dset.num_obs)
    temp = np.empty(dset.num_obs)
    dt = np.empty(dset.num_obs)
    tm = np.empty(dset.num_obs)
    e = np.empty(dset.num_obs)
    ah = np.empty(dset.num_obs)
    aw = np.empty(dset.num_obs)
    mh = np.empty(dset.num_obs)
    mw = np.empty(dset.num_obs)
    la = np.empty(dset.num_obs)
    undu = np.empty(dset.num_obs)

    mjd = dset.time.utc.mjd
    lat, lon, height = dset.site_pos.pos.llh.T
    zd = dset.site_pos.zenith_distance

    # Determine GPT2W values for each observation by interpolating between two unique
    # daily solutions
    for obs in range(dset.num_obs):
        # Start 'gpt2.f' day-by-day in folder where 'gpt2_5.grd' is placed and carry out
        # linear interpolation
        (press[obs], temp[obs], dt[obs], tm[obs], e[obs], ah[obs], aw[obs], la[obs], undu[obs]) = gpt2w_wrapper(
            mjd[obs], [lat[obs]], [lon[obs]], [height[obs]]
        )

        # Determine mapping function values based on coefficients 'ah' and 'aw'
        mh[obs], mw[obs] = iers.vmf1_ht(ah[obs], aw[obs], mjd[obs], lat[obs], height[obs], zd[obs])

    return press, temp, dt, tm, e, mh, mw, la, undu


def gpt2w_wrapper(mjd, lat, lon, hell):
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
        lat (list):           Array with latitude for each station in [rad].
        lon (list):           Array with longitude for each station in [rad].
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
    nstat = len(lat)  # Number of stations
    it = 0  # Use of time variations (annual and semiannual terms)

    # TODO: Case not handled if several stations are included in Dateset. Is that the case for VLBI?

    if not (len(lat) == len(lon) == len(hell)):
        log.fatal("Length of latitude, longitute and ellipsoidal height array is not equal.")

    # Change directory so that gpt2w.f can read the gpt2_5.grd-file in the GPT2w source directory
    current_dir = os.getcwd()
    os.chdir(config.files.path(ext_gpt2w.__name__))

    # Loop over all unique dates (rounded to integer value)
    for date in _rounded_dates(mjd):

        # Check if date is already included in cache
        if date not in _GPT2W:
            _GPT2W[date] = np.array(ext_gpt2w.gpt2_1w(date, lat, lon, hell, nstat, it)).reshape(-1)
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
    #model = config.tech.get("zenith_hydrostatic_delay", section=MODEL, default="saastamoinen").str

    #if model == "saastamoinen":
    # TODO: model needs to in a different configuration parameter than zenith_hydrostatic_delay
    pressure = saastamoinen_pressure(zhd, latitude, height)
    #else:
    #    log.fatal(
    #        f"Zenith troposphere delay definition {model!r} is not correct in configuration file. "
    #        "It should be 'saastamoinen'"
    #    )

    return pressure


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


def site_pressure(dset):
    """Get atmospheric pressure from local site measurements

    If local atmospheric pressure measurements on a site are not available an alternative model given in configuration
    file is used to determine atmospheric pressure.

    TODO:
        So far only gridded VMF1 model is used, if local pressure data are not available.  Which alternative is used,
        should be decided via the configuration file. How to check after an alternative model in configuration file?
        model_list = config.tech[MODEL].meteorological_data.list???

    Args:
        dset (Dataset): A Dataset containing model data.

    Returns:
        numpy.ndarray: Atmospheric pressure for each observation in [hPa]
    """
    pressure = np.zeros(dset.num_obs)

    i_given = np.zeros(dset.num_obs, dtype=bool)
    if "pressure" + (dset.default_field_suffix or "") in dset.fields:
        i_given[np.logical_not(np.isnan(dset.pressure))] = True
        pressure[i_given] = dset.pressure[i_given]

    i_missing = np.logical_not(i_given)
    if i_missing.any():
        log.warn(f"No pressure data for some epochs in dataset. Using VMF1 station data as backup")
        pressure[i_missing] = vmf1_station_pressure(dset)[i_missing]

    return pressure


def vmf1_gridded_pressure(dset):
    """Calulates VMF1 gridded atmospheric pressure

    Equation (9.11) in IERS conventions 2010 :cite:`iers2010` is converted to 'pressure' and the VMF1 gridded
    atmospheric pressure can be determined with help of the gridded zenith hydrostatic delays from VMF1. The gridded
    VMF1 pressure values needs to be rescaled to actual station height. The method is described with Eq. (4) in Kouba
    :cite:`kouba2007`.

    Args:
        dset (Dataset):    Model data.

    Returns:
        numpy.ndarray:  Atmospheric pressure for each observation in [hPa].
    """
    # Get gridded VMF1 data
    vmf1 = apriori.get("vmf1_grid", time=dset.time)
    lat, lon, height = dset.site_pos.pos.llh.T
    try:
        grid_zhd = vmf1["zh"](dset.time, lon, lat)  # Interpolation in time and space in VMF1 grid
        grid_height = vmf1["ell"](lon, lat, grid=False)
        grid_pressure = pressure_zhd(grid_zhd, lat, grid_height)

        # Rescale gridded pressure to station pressure
        pressure = pressure_height_correction(grid_pressure, grid_height, height)
    except KeyError:
        log.warn("No VMF1 grid data available. Unable to compute pressure value. Setting pressure to 1 atm (1013.25hPa). ")
        pressure = np.full(len(lat), fill_value=1013.25)
    return pressure

def vmf1_station_pressure(dset):
    """Calulates atmospheric pressure from station dependend VMF1 files

    Args:
        dset (Dataset):    Model data.

    Returns:
        numpy.ndarray:  Atmospheric pressure for each observation in [hPa].
    """
    vmf1 = apriori.get("vmf1_station", time=dset.time)
    pressure = np.zeros(dset.num_obs)

    for sta in dset.unique("station"):
        idx = dset.station == sta
        try:
            pressure[idx] = vmf1[sta]["pressure"](dset.time.mjd[idx])
        except KeyError:
            # Station is missing, use grid as backup
            log.warn(f"No data for station {sta} in VMF1 station data. Using VMF1 grid as backup")
            pressure_grid = vmf1_gridded_pressure(dset)
            pressure[idx] = pressure_grid[idx]
    return pressure

def vmf1_gridded_mapping_function(dset):
    """Calculates VMF1 hydrostatic and wet mapping functions based on gridded VMF1 files

    This routine determines the VMF1 (Vienna Mapping Functions 1) described in Boehm et al. :cite:`boehm2006a` and
    Kouba :cite:`kouba2007` and uses the 'vmf1_ht.f' from the IERS software library :cite:`iers2010`.

    Args:
        dset (Dataset):    Model data.

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
    vmf1 = apriori.get("vmf1_grid", time=dset.time)

    mh = np.empty(dset.num_obs)
    mw = np.empty(dset.num_obs)

    lat, lon, height = dset.site_pos.pos.llh.T
    zd = dset.site_pos.zenith_distance
    idx_missing = np.zeros(dset.num_obs, dtype=bool)

    for obs in range(dset.num_obs):
        try:
            mh[obs], mw[obs] = iers.vmf1_ht(
                vmf1["ah"](dset.time[obs], lon[obs], lat[obs]),
                vmf1["aw"](dset.time[obs], lon[obs], lat[obs]),
                dset.time.utc.mjd_int[obs],
                lat[obs],
                height[obs],
                zd[obs],
            )
        except KeyError:
            idx_missing[obs] = True
    
    if idx_missing.any():
        log.warn(f"No VMF1 grid data available. Using GMF mapping function.")
        mh_gmf, mw_gmf = gmf_mapping_function(dset)
        mh[idx_missing] = mh_gmf[idx_missing]
        mw[idx_missing] = mw_gmf[idx_missing]

    return mh, mw

def vmf1_station_mapping_function(dset):
    """Calculates VMF1 hydrostatic and wet mapping functions based on station VMF1 files

    This routine determines the VMF1 (Vienna Mapping Functions 1) described in Boehm et al. :cite:`boehm2006a` and
    Kouba :cite:`kouba2007` and uses the 'vmf1_ht.f' from the IERS software library :cite:`iers2010`.

    Args:
        dset (Dataset):    Model data.

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
    vmf1 = apriori.get("vmf1_station", time=dset.time)

    mh = np.empty(dset.num_obs)
    mw = np.empty(dset.num_obs)

    lat = dset.site_pos.pos.llh.lat
    zd = dset.site_pos.zenith_distance
    sta = dset.station
    mjd = dset.time.mjd
    missing_idx = np.zeros(dset.num_obs, dtype=bool)

    for obs in range(dset.num_obs):
        try:
            mh[obs], mw[obs] = iers.vmf1(
                vmf1[sta[obs]]["ah"](mjd[obs]),
                vmf1[sta[obs]]["aw"](mjd[obs]),
                mjd[obs],
                lat[obs],
                zd[obs],
            )
        except KeyError:
            # Station is missing, use grid as backup
            missing_idx[obs] = True

    if missing_idx.any():
        log.warn(f"No data for some epochs in VMF1 station data. Using VMF1 grid as backup")
        mh_grid, mw_grid = vmf1_gridded_mapping_function(dset)
        mh[missing_idx] = mh_grid[missing_idx]
        mw[missing_idx] = mw_grid[missing_idx]

    return mh, mw


def vmf1_gridded_zenith_wet_delay(dset):
    """Calculates zenith wet delay based on gridded zenith wet delays from VMF1

    Uses gridded zenith wet delays from VMF1, which are rescaled from the gridded height to actual station height by
    using Equation(5) described in Kouba :cite:`kouba2007`.

    Args:
        dset (Dataset):    Model data.

    Returns:
        numpy.ndarray:     Zenith wet delay for each observation in [m]
    """
    # Get gridded VMF1 data
    vmf1 = apriori.get("vmf1_grid", time=dset.time)

    lat, lon, height = dset.site_pos.pos.llh.T
    try:
        grid_zwd = vmf1["zw"](dset.time, lon, lat)  # Interpolation in time and space in VMF1 grid
        grid_height = vmf1["ell"](lon, lat, grid=False)

        # Zenith Wet delay. Eq. (5) in Kouba :cite:`kouba2007`
        zwd = grid_zwd * np.exp(-(height - grid_height) / 2000)
    except KeyError:
        log.warn(f"No VMF1 grid data available. Setting zenith wet delay to 0.")
        zwd = np.zeros(dset.num_obs)

    return zwd

def vmf1_station_zenith_wet_delay(dset):
    """Calculates zenith wet delay based on station zenith wet delays from VMF1

    Args:
        dset (Dataset):    Model data.

    Returns:
        numpy.ndarray:     Zenith wet delay for each observation in [m]
    """
    vmf1 = apriori.get("vmf1_station", time=dset.time)
    zwd = np.zeros(dset.num_obs)

    for sta in dset.unique("station"):
        idx = dset.station == sta
        try:
            zwd[idx] = vmf1[sta]["zw"](dset.time.mjd[idx])
        except KeyError:
            # Station is missing, use grid as backup
            log.warn(f"No data for station {sta} in VMF1 station data. Using VMF1 grid as backup")
            zwd_grid = vmf1_gridded_zenith_wet_delay(dset)
            zwd[idx] = zwd_grid[idx]

    return zwd

def vmf1_gridded_zenith_hydrostatic_delay(dset):
    """Calculates zenith hydrostatic delay based on gridded zenith wet delays from VMF1

    Uses gridded zenith wet delays from VMF1, which are rescaled from the gridded height to actual station height by
    using Equation (3) and (4) described in Kouba :cite:`kouba2007`.

    Args:
        dset (Dataset):    Model data.

    Returns:
        numpy.ndarray:     Zenith hydrostatic delay for each observation in [m]
    """
    # Gridded pressure rescaled to station height
    pressure = vmf1_gridded_pressure(dset)
    # Zenith hydrostatic delay. Eq. (4) in Kouba :cite:`kouba2007`
    lat, _, height = dset.site_pos.pos.llh.T
    zhd = saastamoinen_zenith_hydrostatic_delay(pressure, lat, height)

    return zhd

def vmf1_station_zenith_hydrostatic_delay(dset):
    """Calculates zenith hydrostatic delay based on station zenith hydrostatic delays from VMF1

    Args:
        dset (Dataset):    Model data.

    Returns:
        numpy.ndarray:     Zenith hydrostatic delay for each observation in [m]
    """
    vmf1 = apriori.get("vmf1_station", time=dset.time)
    zhd = np.zeros(dset.num_obs)

    for sta in dset.unique("station"):
        idx = dset.station == sta
        try:
            zhd[idx] = vmf1[sta]["zh"](dset.time.mjd[idx])
        except KeyError:
            # Station is missing, use grid as backup
            log.warn(f"No data for station {sta} in VMF1 station data. Using VMF1 grid as backup")
            zhd_grid = vmf1_gridded_zenith_hydrostatic_delay(dset)
            zhd[idx] = zhd_grid[idx]

    return zhd


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
