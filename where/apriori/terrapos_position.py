"""Get Terrapos position Dataset by reading Terrapos position output file

Example:
from where import apriori
dset = apriori.get('terrapos_position') # File key 'terrapos_output_position' has to be defined in configuration file
                                        # files.conf.

# Rename the dataset and write it to file.
dset.write_as(tech=tech, stage=stage, dataset_name=dataset_name, dataset_id=dataset_id)

Description:
------------
Read Terrapos position output file and enable correct time and position format to create a Dataset.



"""
# External library imports
import numpy as np

# Where imports
from where.lib import gnss
from where.lib import plugins
from where.lib.time import Time
from where import parsers


@plugins.register
def get_terrapos_position():
    """Get Terrapos position Dataset by reading Terrapos position output file

    Returns:
        dset (Dataset): Dataset with following fields:

        ====================  ==================  =================================================================
         Field                 Type                Description
        ====================  ==================  =================================================================
         gpsweek               numpy.ndarray       GPS week
         gpssec                numpy.ndarray       Seconds of GPS week
         head                  numpy.ndarray       Head in [deg]
         height                numpy.ndarray       Ellipsoidal height in [m]
         lat                   numpy.ndarray       Latitude in [deg]
         lon                   numpy.ndarray       Longitude in [deg]
         num_sat               numpy.ndarray       Number of satellites
         pdop                  numpy.ndarray       Position Dilution of Precision (PDOP)
         pitch                 numpy.ndarray       Pitch in [deg]
         reliability_east      numpy.ndarray       East position external reliability in [m] #TODO: Is that correct?
         reliability_height    numpy.ndarray       Height position external reliability in [m] #TODO: Is that correct?
         reliability_north     numpy.ndarray       North position external reliability in [m] #TODO: Is that correct?
         roll                  numpy.ndarray       Roll in [deg]
         sigma_east            numpy.ndarray       Standard deviation of East position in [m] #TODO: Is that correct?
         sigma_height          numpy.ndarray       Standard deviation of Height position in [m] #TODO: Is that correct?
         sigma_north           numpy.ndarray       Standard deviation of North position in [m] #TODO: Is that correct?
         site_pos              PositionTable       PositionTable object with given station coordinates
         time                  TimeTable           Observation time given as TimeTable object
        ====================  ==================  =================================================================
    """
    dset = parsers.parse_key("terrapos_output_position").as_dataset()
    dset.add_time("time", val=_get_time(dset), scale="gps")
    dset.add_position("site_pos", time="time", itrs=_get_site_pos(dset))

    return dset


def _get_site_pos(dset):
    """Determine site position by converting given longitude, latitude and height for a station to geocentric
       cartesian coordinates
    """
    # TODO hjegei: Workaround -> better would it be if Position object can handle LLH as input format!!!
    x, y, z = gnss.llh2xyz(np.deg2rad(dset.lat), np.deg2rad(dset.lon), dset.height)
    return np.stack((x, y, z), axis=1)


def _get_time(dset):
    """Determine time field
    """
    # TODO hjegei: Workaround -> better would it be if Time object can handle gpsweek as input format!!!
    jd_day, jd_frac = gnss.gpssec2jd(dset.gpsweek, dset.gpssec)
    return Time(val=jd_day, val2=jd_frac, format="jd", scale="gps")
