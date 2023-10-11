"""a GNSS High Accuracy Service (Galileo) pipeline for generation of SP3 files

Description:
------------

SP3 files are generated, whereby broadcast navigation messages are read and corrected by corresponding HAS messages.
The interval in the SP3 files depends on simulated observation epochs. For these observation epochs the orbit and
clock are determined. 

Following steps are normally carried out:
    
    1. Setup stage
       Generate simulated dataset, which includes the fields 'time', 'system' and 'satellite' for a given sampling
       rate.
    2. Edit stage
       Data are edited in dependency on defined removers (_ignore_satellites) in configuration file.
    3. Orbit stage
        a.) Read Galileo HAS orbit corrections and remove observations, which should be used due to
            unavailable HAS messages.
        b.) Remove observation epochs which are below HAS receiver reception time or exceeding validity
            length of Galileo HAS orbit correction messages.
        c.) Read broadcast navigation messages and remove duplicated messages, observation epochs exceeding validity
            length of broadcast navigation message record and unhealthy satellites.
        d.) Keep only observations, where IOD does match between HAS and broadcast navigation message.
        e.) Apply antenna phase center offset in case SP3 file should be refered to Center of Mass.
        f.) Apply Galileo orbit correction on broadcast ephemeris.
    
        g.) Read Galileo HAS clock corrections and remove observations, which should be used due to
            unavailable HAS messages.
        h.) Remove observation epochs which are below HAS receiver reception time or exceeding validity
            length of Galileo HAS clock correction messages.
        i.) Apply Galileo orbit correction on broadcast ephemeris.
    4. Write stage
       Write output files in dependency on defined writers.

"""

# Standard library imports
from datetime import datetime, timedelta

# External library imports
from typing import Dict, List, Tuple, Union

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where import cleaners
from where import writers
from where.lib import config
from where.lib import log
from where.lib import util

# The name of this technique
PIPELINE = __name__.split(".")[-1]

TEST = False


@plugins.register_named("options")
def options() -> Tuple[str]:
    """Command line options that can be used to specify this technique

    Returns:
        Strings specifying command line options.
    """
    return "--gnss_has_sp3",


@plugins.register_named("get_args")
def get_args(rundate: datetime.date, input_args: Union[None, List[str]]=None) -> List[str]:
    """Convert where_runner arguments to where arguments for given date

    Args:
        rundate:    The model run date.
        input_args: Input arguments.

    Returns:
        Strings with names of available sessions.
    """
    return [" ".join(input_args)] if input_args else list()


@plugins.register_named("file_vars")
def file_vars(file_vars=None) -> Dict[str, str]:
    """File variables that will be available during the running of this technique

    In addition, date and analysis variables are available.

    Returns:
        File variables special for this technique.
    """
    station = config.tech.station.str

    # Determine correct interval based on given sampling rate
    sampling_rate = config.tech.sampling_rate.int
    if sampling_rate < 59: # seconds
        interval = str(sampling_rate).zfill(2) + "S"
    elif sampling_rate > 59 & sampling_rate < 3599: # seconds
        interval = str(int(sampling_rate/60)).zfill(2) + "M"
    elif sampling_rate > 3599 & sampling_rate < 86399: # seconds
        interval = str(int(sampling_rate/3600)).zfill(2) + "H"
    elif sampling_rate > 86399 & sampling_rate < 8639999: # seconds
        interval = str(int(sampling_rate/3600/24)).zfill(2) + "D"
    elif sampling_rate > 8639999:
        log.fatal(f"Interval based sampling_rate {sampling_rate} is not defined.")
        
    return dict(interval=interval, station=station.lower(), STATION=station.upper())


#
# SETUP SIMULATED OBSERVATION EPOCHS
#
@plugins.register
def setup(stage: str, dset: "Dataset") -> None:
    """Setup Where dataset for getting observation epochs

    The dataset includes the fields 'time', 'system' and 'satellite' like:
        <epoch>  <system>  <satellite>

    Example:

    EPOCH                                  SYSTEM   SATELLITE
    ---------------------------------------------------------------
    datetime.datetime(2016, 3, 1, 1, 0)     'G'      'G01'
    datetime.datetime(2016, 3, 1, 1, 0)     'G'      'G02'
    ...
    datetime.datetime(2016, 3, 1, 1, 0)     'E'      'E01'
    datetime.datetime(2016, 3, 1, 1, 0)     'E'      'E02'
    ...
    datetime.datetime(2016, 3, 1, 1, 5)     'G'      'G01'
    datetime.datetime(2016, 3, 1, 1, 5)     'G'      'G02'
    ...

    The time stamps used for the observation epochs are based on the sampling rate. The 'rundate' at 00:00 o'clock in 
    GPS time scale defines the start time. The sampling rate given in the configuration file defines the time offset 
    added to start time until 24:00 o'clock in GPS time scale is reached.

    The GNSS satellites are defined in the configuration file, which should be used for generation of the SP3 file.

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    # Initialize variables
    time = list()
    used_satellites = list()
    used_systems = list()
    dset_satellites = list()
    dset_systems = list()

    dset.vars.update(file_vars())

    # Get options from configuration file
    sampling_rate = config.tech.sampling_rate.float
    systems = config.tech.systems.list
    satellites = config.tech.satellites.list
    dset.meta["sampling_rate"] = sampling_rate

    # Select satellites for further orbit comparison
    for satellite in satellites:
        if satellite.startswith(tuple(systems)):
            used_satellites.append(satellite)
            used_systems.append(satellite[0])  # Needed for generation of Dataset

    # Prepare generation of Dataset
    hour = 0
    minute = 0
    second = 0
    start_time = datetime(
        dset.analysis["rundate"].year,
        dset.analysis["rundate"].month,
        dset.analysis["rundate"].day,
        hour,
        minute,
        second,
    )
    time_to_save = start_time
    num_satellite = len(used_satellites)
    while time_to_save < start_time + timedelta(days=1):
        time.extend([time_to_save] * num_satellite)
        time_to_save += timedelta(seconds=sampling_rate)
        dset_satellites.extend(used_satellites)  # All satellites has to be given in each epoch
        dset_systems.extend(used_systems)  # All systems identifiers has to be given in each epoch

    # Generate Dataset
    dset.num_obs = len(time)
    dset.add_time("time", val=time, scale="gps", fmt="datetime", write_level="operational")
    dset.add_text("satellite", val=dset_satellites, write_level="operational")
    dset.add_text("system", val=dset_systems, write_level="operational")

    # Write Dataset to file
    if util.check_write_level("analysis"):
        dset.write_as(stage="setup")


#
# EDIT DATA
#
@plugins.register
def edit(stage: str, dset: "Dataset") -> None:
    """Edit the Galileo HAS data.

    Args:
        stage:      Name of current stage.
        dset:       A dataset containing the data.
    """
    cleaners.apply_removers("removers", dset)
    
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)
    
    
#
# ORBIT DETERMINATION
#
@plugins.register
def orbit(stage: str, dset: "Dataset") -> None:
    """Determine broadcast orbit and clocks corrected with HAS corrections

    Args:
        stage:  Name of current stage.
        dset:   A dataset containing the data.
    """
    satellite_origin = config.tech.satellite_origin.str
      
    # Add HAS orbit correction. 
    #
    # Note: This step has to be done first to match the broadcast navigation message to the HAS messages and to add the 
    #       "sat_posvel", "delay.gnss_satellite_clock" and "delay.gnss_relativistic_clock" fields to dataset.
    orbit_has = apriori.get(
        "orbit", 
        rundate=dset.analysis["rundate"], 
        file_key="gnss_has_orb",
        day_offset=0, 
        apriori_orbit="has",
    ) 
    orbit_has.add_orbit_to_dataset(dset)

    # Add HAS clock correction
    clock_has = apriori.get(
        "orbit", 
        rundate=dset.analysis["rundate"], 
        file_key="gnss_has_clk",
        day_offset=0, 
        apriori_orbit="has",
    )
    clock_has.add_clock_to_dataset(dset)
    
    # Apply HAS orbit and clock correction
    if satellite_origin == "com":
        ant_brdc = apriori.get("gnss_antenna_correction", file_key="antex_brdc") #TODO: Should we IGS ANTEX file? apriori.get("gnss_antenna_correction") 
        
        # Definition of GNSS frequencies for accessing correct PCOs for conversion from APC to CoM:
        #   for broadcast orbits: PCOs of broadcast orbits are based for Galileo on the average of E1/E5/E6 PCOs defined
        #                         by the European GNSS Service Centre and for GPS it is defined by the National Geospatial
        #                         Intelligence Agency (NGA). An ANTEX file is created with these PCOs values, whereby the
        #                         Galileo E1 and respectively GPS L1 is used.
        brdc_sys_freq = {"E": "E1", "G": "L1"} #TODO: Do we need something like that? {"E": "E1_E5b", "G": "L1_L2"} 
        
        pco_sat = ant_brdc.satellite_phase_center_offset(dset, brdc_sys_freq)
        dset.sat_posvel[:] = dset.sat_posvel - pco_sat + dset.has_orbit_correction 
        
    elif satellite_origin == "apc":
        dset.sat_posvel[:] = dset.sat_posvel + dset.has_orbit_correction
    else:
        log.fatal(f"Satellite origin {satellite_origin} is unknown, either CoM (Center of Mass) or APC (Antenna " 
                  f"phase center) can be choosen.")
        
    dset.delay.gnss_satellite_clock[:] = dset.delay.gnss_satellite_clock - dset.has_clock_correction
                   
    if util.check_write_level("analysis"):
        dset.write_as(stage=stage)    
        
    
#
# WRITE RESULTS
#
@plugins.register
def write(stage: str, dset: "Dataset") -> None:
    """Write results to file.

    Write results to file. This uses the writers framework which calls different writers depending on the output-field
    in the config-file.

    Args:
        stage (str):          Name of current stage.
        dset (Dataset):       A dataset containing the data.
    """
    writers.write(default_dset=dset)
    if util.check_write_level("operational"):
        dset.write_as(stage="write")
    
    

