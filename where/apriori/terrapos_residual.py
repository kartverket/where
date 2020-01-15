"""Get Terrapos residual Dataset by reading Terrapos residual output file

Example:
from where import apriori
dset = apriori.get('terrapos_residual') # File key 'terrapos_output_residual' has to be defined in configuration file
                                        # files.conf.

# Rename the dataset and write it to file.
dset.write_as(tech=tech, stage=stage, dataset_name=dataset_name, dataset_id=dataset_id)

Description:
------------
Read Terrapos residual output file and enable correct time to create a Dataset.



"""
# Midgard imports
from midgard.dev import plugins

# Where imports
from where.lib import gnss
from where.data.time import Time
from where import parsers


@plugins.register
def get_terrapos_residual():
    """Get Terrapos residual Dataset by reading Terrapos residual output file

    Returns:
        dset (Dataset): Dataset with following fields:

        ====================  ==================  =================================================================
         Field                 Type                Description
        ====================  ==================  =================================================================
         azimuth               numpy.ndarray       Azimuth of satellites in [deg]
         elevation             numpy.ndarray       Elevation of satellites in [deg]
         gpsweek               numpy.ndarray       GPS week
         gpssec                numpy.ndarray       Seconds of GPS week
         residual_code         numpy.ndarray       Code (pseudorange) residuals in [m]
         residual_doppler      numpy.ndarray       Doppler residuals in [m]
         residual_phase        numpy.ndarray       Carrier-phase residuals in [m]
         satellite             numpy.ndarray       Satellite PRN number together with GNSS identifier (e.g. G07)
         system                numpy.ndarray       GNSS identifier
         time                  TimeTable           Observation time given as TimeTable object
        ====================  ==================  =================================================================

    """
    dset = parsers.parse_key("terrapos_output_residual").as_dataset()
    dset.add_time("time", val=_get_time(dset), scale="gps")

    return dset


def _get_time(dset):
    """Determine time field
    """
    # TODO hjegei: Workaround -> better would it be if Time object can handle gpsweek as input format!!!
    jd_day, jd_frac = gnss.gpssec2jd(dset.gpsweek, dset.gpssec)
    return Time(val=jd_day, val2=jd_frac, fmt="jd", scale="gps")
