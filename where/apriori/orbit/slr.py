"""Handling of external SLR orbits

Description:
------------
The module includes a class for handling apriori SLR orbits.

Example:
    from where import apriori

    # Get slr ephemeris object
    orbit = apriori.get('orbit', rundate=rundate, time=time, satellite=satellite, system=system, station=station,
                          apriori_orbit='slr')

    # Write calculated Dataset to file
    orbit.dset.write()


$Revision: 15011 $
$Date: 2018-05-04 16:19:35 +0200 (Fri, 04 May 2018) $
$LastChangedBy: hjegei $

"""
# External imports
from datetime import timedelta

# Where imports
from where.apriori import orbit
from where import data
from where.lib import config
from where.lib import files
from where.lib import log
from where.lib import plugins
from where import parsers


@plugins.register
class Slr(orbit.AprioriOrbit):
    """A class for representing apriori slr orbits

    SP3 orbit files can be read. 

    Attributes:
        day_offset (int):       Day offset used to calculate the day to read.
        dset_orbit (Dataset):   Dataset object, which includes orbits read from SP3 file 
        dset_raw (Dataset):     Dataset object, which includes observations etc.
        file_key (str):         Key to the orbit file defined in files.conf file.
        file_path (pathlib.PosixPath):  File path to SP3 orbit file.
        satellite (tuple):      Satellite number 
        time (Time):            epochs


    Methods:
        _read():                Read orbit file data and save it in a Dataset

    """

    def __init__(self, rundate, time, satellite, system=None, file_key=None, file_path=None, day_offset=6, **kwargs):
        """Set up a new PreciseOrbit object, does not parse any data

        TODO: Remove dependency on rundate, use time to read correct files. (What to do with dataset?)

        Args:
            rundate (date):     Date of model run.
            time (Time):        Time epochs at the satellite for which to calculate the apriori orbit.
            satellite (list):   Strings with names of satellites.
            system (list):      Strings with codes of systems (G, E, R, etc.).
            file_key (str):     Key to the precise orbit file that will be read.
            file_path (pathlib.PosixPath):  File path to SP3 orbit file.
            day_offset (int):   Day offset used to calculate the number of days to read.
        """
        super().__init__(rundate=rundate, time=time, satellite=satellite, system=system)
        self.file_key = file_key
        self.day_offset = day_offset
        self.file_path = file_path

    def _read(self, dset_raw, provider, version):
        """Read SP3 orbit file data and save it in a Dataset

        Naming convention correspond to end of arc, at midnight, hence we add day_offset, 
        which is usually arc_length - 1

        Args:
            dset_raw (Dataset):   Dataset representing raw data from apriori orbit files
            provider:             Str: Orbit provider
            version:              Str: Orbit version
        """
        rundate = dset_raw.rundate
        date_to_read = rundate + timedelta(days=self.day_offset)
        file_vars = config.date_vars(date_to_read)
        file_vars["provider"] = provider
        file_vars["version"] = version

        if self.file_path is None:
            file_path = files.path(self.file_key, file_vars)
        else:
            file_path = self.file_path

        log.debug("Parse precise orbit file {}.", file_path)

        # Generate temporary Dataset with orbit file data
        dset_orbit = data.Dataset(
            rundate=date_to_read, tech=dset_raw.vars["tech"], stage="orbit", dataset_name="", dataset_id=0, empty=True
        )
        parser = parsers.parse(parser_name="orbit_sp3", file_path=file_path, rundate=date_to_read)
        parser.write_to_dataset(dset_orbit)

        dset_orbit.add_to_meta("parser", "file_path", file_path)
        return dset_orbit
