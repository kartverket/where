"""A parser for reading files generated of JRC Galileo HAS decoder

Example:
--------

    from where import parsers
    parser = parsers.parse('gnss_has_decoder', file_path='SEPT147.21_has_cb.csv')

Description:
------------

Reads data from GNSS receiver type file.


"""
# Standard library imports
from pathlib import Path
from typing import Callable, List

# Third party imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.data.time import Time
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.parsers.csv_ import CsvParser

# Where imports
from where.data import dataset3 as dataset
from where.lib import log


@plugins.register
class GnssHasDecoderParser(CsvParser):
    """A parser class for handling of reading files generated of JRC Galileo HAS decoder
    
    The HAS decoder data header line is used to define the keys of the **data** dictionary. The values of the 
    **data** dictionary are represented by the HAS decoder colum values.

    Depending on the HAS decoder file following dataset fields can be available:
    
    | Field                    | Type   | Tbl | Description                                                           |
    |--------------------------|--------------|-----------------------------------------------------------------------|
    | av_flag                  | CB, CP |   7 | Availablity flag related to code or phase bias, which represents      |
    |                          |        |     | number of available signals. The number of available signals          |
    |                          |        |     | corresponds to the number of given biases. If av_flag = 0, then no    |
    |                          |        |     | data are available.                                                   |
    | code_bias                | CB     |   7 | Code bias of satellite in [m]                                         |
    | delta_clock_c0           | CLK    |   7 | Delta clock C0 correction for satellite in [m].                       |
    | delta_cross_track        | ORB    |   7 | Delta cross-track correction for satellite in [m]. Not a number (NaN),|
    |                          |        |     | if data are not available.                                            |
    | delta_in_track           | ORB    |   7 | Delta in-track correction for satellite in [m]. Not a number (NaN),   |
    |                          |        |     | if data are not available.                                            |
    | delta_radial             | ORB    |   7 | Delta radial correction for satellite in [m]. Not a number (NaN),     |
    |                          |        |     | if data are not available.                                            |
    | gnssID                   | ALL    |   8 | GNSS identifier (0: GPS, 2: Galileo)                                  |
    | gnssIOD                  | ALL    |   7 | Reference IOD for which orbital corrections are determined. For CB,   |
    |                          |        |     | CP and CLK HAS messages is gnssIOD derived from gnssIOD information   |
    |                          |        |     | in ORB files by using IOD, which connects HAS messages to the same    |
    |                          |        |     | IOD set.                                                              |
    | IOD                      | ALL    |   6 | Issue of data extracted from HAS header                               |
    | multiplier               | CLK    |   7 | Multiplier for all delta clock C0 corrections of GNSS ID              |
    | phase_bias               | CP     |   7 | Phase bias of satellite in [cycle]                                    |
    | phase_discontinuity_ind  | CP     |   7 | Phase discontinuity counter                                           |
    | PRN                      | ALL    |   9 | Galileo SVID or GPS PRN                                               |
    | signal                   | CB, CP |  10 | Signal identifier                                                     |
    | status                   | CLK    |   7 | Delta clock C0 correction status:                                     |
    |                          |        |     |    0 - data ok                                                        |
    |                          |        |     |    1 - data not available                                             |
    |                          |        |     |    2 - satellite shall not be used                                    |
    | ToH                      | ALL    |   6 | Time of Hour of the message information related to Galileo System     |
    |                          |        |     | Time (GST) in [s]  extracted from HAS header                          |
    | ToW                      | ALL    |     | Time of Week in [s] extracted from navigation message                 |
    | validity                 | ALL    |  12 | Validity interval in [s] starting with ToH                            |
    | WN                       | ALL    |     | GPS week number extrated from navigation message                      |
    
    Following time fields are added in addition:

    | Field                    |  Type  | Description                                                                 |
    |--------------------------|--------------------------------------------------------------------------------------|
    | time                     | Time   | Receiver reception time of HAS message                                      |
    | tom                      | Time   | Reference time of HAS message                                               |

    """
    
    #
    # SETUP POSTPROCESSORS
    #
    def setup_postprocessors(self) -> List[Callable[[], None]]:
        """List steps necessary for postprocessing
        """
        return [
            #self._clean_messages, #Note: Postprocessor has to be carried out before _add_time().
            self._add_time,
        ]
 
    
    def _add_time(self) -> None:
        """Add receiver reception time of HAS messages and reference time of HAS message to data
        
        #TODO: HAS format should be improved, that also GPS week is included in HAS message file.
        """
        num_obs = len(self.data["ToW"])
       
        # Determine receiver reception time of HAS message
        time = Time(
            val=self.data["WN"], 
            val2=self.data["ToW"], 
            scale="gps", 
            fmt="gps_ws",
        )

        # Determine time of HAS message (TOM), which is the reference time of the HAS message
        tom = np.zeros(num_obs)
        tom_hour = np.floor(time.gps_ws.seconds /3600) * 3600  # Hour of time of message in GPS seconds
        
        # Handling of ToH in HAS message related to Eq. 24 of Galileo HAS SIS ICD v1.5
        idx = tom_hour + self.data["ToH"] <= time.gps_ws.seconds
        
        if np.any(idx):
            tom[idx] = tom_hour[idx] + self.data["ToH"][idx]
            
        if np.any(~idx):  # ~idx is inversion of numpy boolean array idx 
            tom[~idx] = (tom_hour[~idx] - 3600) + self.data["ToH"][~idx]
                      
        tom = Time(
            val=self.data["WN"],
            val2=tom,
            scale="gps", 
            fmt="gps_ws",
        )
            
        # Add receicer reception time and time of HAS message (reference time)
        self.data["time"] = time
        self.data["tom"] = tom
        
        
    def _clean_messages(self) -> None:
        """Clean HAS messages
        
        Remove HAS messages, which include failures. At the moment are only HAS clock correction handled.
        """
        #TODO: Check decoded HAS message data on day 22.03.2022 (doy: 081)
        if "multiplier" in self.data.keys():
            df = pd.DataFrame.from_dict(self.data)

  
    def as_dataset(self) -> "Dataset":
            """Return the parsed data as a Dataset
    
            Returns:
                A dataset containing the data.
            """
            # GNSS signal definition of HAS message based on Table 10 for Galileo and GPS
            signal_def_gal = {
                0: "E1-B",
                1: "E1-C",
                2: "E1-B+E1-C",
                3: "E5a-I",
                4: "E5a-Q",
                5: "E5a-I+E5a-Q",
                6: "E5b-I",
                7: "E5b-Q",
                8: "E5b-I+E5b-Q", 
                9: "E5-I",
               10: "E5-Q",
               11: "E5-I+E5-Q",
               12: "E6-B",
               13: "E6-C",
               14: "E6-B+E6-C",
            }
            
            signal_def_gps = {
                0: "L1 C/A",
                3: "L1C(D)",
                4: "L1C(P)",
                5: "L1C(D+P)",
                6: "L2 CM",
                7: "L2 CL",
                8: "L2 CM+CL",
                9: "L2 P",
               11: "L5 I", 
               12: "L5 Q",
               13: "L5 I+L5 Q",
            }
            
            # GNSS definition of HAS message based on Table 8
            system_def = {
                0: "G", # GPS
                2: "E", # Galileo
            }
            
            # Define unit of float fields
            unit_def = {
                "code_bias": "m",
                "delta_clock_c0": "m",
                "delta_cross_track": "m",
                "delta_in_track": "m",
                "delta_radial": "m",
                "phase_bias": "", #TODO: cycle!!!
                "validity": "s",
            }
                        
            # Initialize dataset
            dset = dataset.Dataset()
            if not self.data:
                log.warn("No data in {self.file_path}.")
                return dset
            dset.num_obs = len(self.data["ToW"])
            
            # Add time fields
            for field in ["time", "tom"]:
                dset.add_time(
                    name=field,
                    val=self.data[field],
                    write_level="operational",
                )
    
            # Add system field based on gnssID column
            dset.add_text("system", val=[system_def[value] for value in self.data["gnssID"]])
    
            # Add satellite field based on PRN column
            prn_data = []
            for idx, prn in enumerate(self.data["PRN"]):
                prn_data.append(system_def[self.data["gnssID"][idx]] + str(prn).zfill(2))
    
            dset.add_text("satellite", val=prn_data)
    
            # Add signal field based on signal column
            if "signal" in self.data.keys():
                values = np.empty(dset.num_obs, dtype="object_")
                for sys in dset.unique("system"):
                    idx = dset.filter(system=sys)
                    if sys == "E":
                        values[idx] = [signal_def_gal[v] for v in self.data["signal"][idx]]
                    elif sys == "G":
                        values[idx] = [signal_def_gps[v] for v in self.data["signal"][idx]]
                    else:
                        log.fatal(f"GNSS '{sys}' is not defined for HAS messages. Only 'E=2'-Galileo and 'G=0'-GPS " 
                                  f"are defined.")
                dset.add_text("signal", val=values)
                
            # Add age of data
            dset.add_float(
                        "age_of_data",
                        val = ((dset.time.gps.mjd - dset.tom.gps.mjd) * Unit.day2second), 
                        unit = "second",
            ),
    
            # Define fields to save in dataset
            remove_fields = {"gnssID", "PRN", "signal", "time", "ToH", "tom", "ToW", "WN"}
            fields = set(self.data.keys()) - remove_fields
    
            # Add float fields
            for field in fields:
                unit = unit_def[field] if field in unit_def else ""
                dset.add_float(field.lower(), val=self.data[field], unit=unit)
                
            # Add information about file type in meta variable
            file_type_def = {
                "cb": "code bias",
                "clk": "clock",
                "cp": "phase bias",
                "orb": "orbit",
            }
            file_type = Path(self.meta["__data_path__"]).stem.split('_')[-1]
            dset.meta["file_type"] = file_type_def[file_type] if file_type in file_type_def else "unknown"
    
            return dset
