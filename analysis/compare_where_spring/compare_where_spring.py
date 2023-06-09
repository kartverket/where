#!/usr/bin/env python3
"""Read Spring CSV output files and convert them to Where dataset

Description:
------------

"""

# Standard library imports
from datetime import date
import pathlib
import subprocess
from typing import List

# External library imports
import numpy as np

# Midgard imports
from midgard.data import position
from midgard.data._time import GpsTime, UtcTime

# Where imports
from where.data import dataset3 as dataset
from where import parsers
from where.lib import config
from where.lib import log


def _concatenate_fields(dset: "Dataset", fields: List[str]) -> List[str]:
    """Concatenate fields to string lines

    Args:
        dset:    A dataset containing data.
        fields:  Name of fields to be concatenated

    TODO: Write a general routine like Dataframe function "to_string()".

    Returns:
        List with string lines, whereby each line represents concatenated fields
    """
    concatenated_fields = []
    for values in dset.values(*fields):
        line = ""
        for value in values:
            if isinstance(value, (GpsTime, UtcTime)):  # TODO: Check if more time datatypes should be defined.
                line = f"{line}{value.gps.isot}"
            else:
                line = f"{line}{value}"

        concatenated_fields.append(line.strip())

    return concatenated_fields


def _generate_where_dset(rundate, pipeline, station, stage, system):
    """Generate Where dataset

    Start Where analysis with Where and read Dataset.
    """
    program = "where"

    # Start processing with Where
    where_args = [
        str(rundate.year),
        str(rundate.month),
        str(rundate.day),
        "--" + pipeline,
        "--station=" + station,
        "--systems=" + system,
        "--sampling_rate=30", #MURKS =30
        "--estimate_epochwise=True",
        "--elevation:cut_off=5",
        "--elevation_weighting=none",
        "--brdc_block_nearest_to=transmission_time:positive",
        "--site=",
        "-N",
        "-D",
        "-T",
    ]
    
    if pipeline == "gnss":
        
        if system == "G":
            where_args.extend([
                  #"--gnss_ionosphere:model=klobuchar",
                  #"--profile=cnes_lnav_l1",  
                  "--profile=cnes_lnav_l1l2",  
            ])
        elif system == "E":
            where_args.extend([
                  "--gnss_ionosphere:model=nequick",
                  # "--gnss_ionosphere:model=klobuchar",
                  "--profile=cnes_inav_e1",
                  #"--profile=cnes_inav_e1e5b",
                  #"--profile=cnes_fnav_e1e5a",
            ])

    elif pipeline == "gnss_vel":
        where_args.extend(["--gnss_select_obs:obs_code=code,doppler"])
        
    log.info("Start {}: {} {}".format(pipeline.upper(), program, " ".join(where_args)))
    #process = subprocess.run([program] + where_args)
    #if process.returncode:
    #    log.error("Where {} failed with error code {} ({})", pipeline.upper(), process.returncode, " ".join(process.args))

    # Read and write GNSS Dataset
    dset = dataset.Dataset.read(rundate=rundate, pipeline=pipeline, stage=stage, station=station, label="None", id="")

    # Add sat_pos field
    dset.add_position("sat_pos", time=dset.time, system="trs", val=dset.sat_posvel.trs.pos)
    dset.write_as(pipeline=pipeline, stage="where")

    return dset


def _merge_diff_dset(dset_diff_merge, dset_diff):

    merge_by = ["time"]

    for field in dset_diff.fields:
        if field not in dset_diff_merge.fields:

            table = dset_diff._data[dset_diff._fields[field]].__module__.split(".")[-1].split("_")[0]

            # Get common dataset data
            hash_merge = _concatenate_fields(dset_diff_merge, merge_by)
            hash_diff = _concatenate_fields(dset_diff, merge_by)
            hash_both = sorted(set(hash_merge) & set(hash_diff))
            idx_merge = [hash_merge.index(h) for h in hash_both]

            # Decimate merged dataset to correct number of observations
            dset_diff_merge.subset(idx_merge)

            # Add field
            if table == "float":
                shape = (3,) if len(dset_diff[field].shape) == 2 else (1,)
                dset_diff_merge.add_float(field, val=dset_diff[field], shape=shape)

            elif table == "text":
                dset_diff_merge.add_text(field, val=dset_diff[field])

            elif table == "time":
                dset_diff_merge.add_time(field, val=dset_diff[field], scale=dset_diff[field].scale)

            elif table == "position" or table == "posvel":
                getattr(dset_diff, f"add_{table}")(field, itrs=dset_diff[field].itrs)


def _read_spring_files(file_, rundate, pipeline, station):

    stage = "spring" + str(file_.stem)[0:3].lower()

    # Read already existing dataset
    file_vars = dict(
                stage=stage,
                station=station,  
                label="None", 
                id=""
    )
    if config.files.path("dataset", file_vars=file_vars).exists():

        dset = dataset.Dataset.read(rundate=rundate, pipeline=pipeline, **file_vars)

        return dset

    # Generate dataset
    if file_.exists():

        # Read file
        p = parsers.parse_file(parser_name="spring_csv", file_path=file_)
        dset = p.as_dataset()

    else:
        log.warn(f"{file_} does not exists.")
        return

    # Add HPE and VPE to dataset
    if (
        "site_pos_vs_ref_east" in dset.fields
        and "site_pos_vs_ref_north" in dset.fields
        and "site_pos_vs_ref_vertical" in dset.fields
    ):
        dset.add_float(
            "hpe", val=np.sqrt(dset.site_pos_vs_ref_east ** 2 + dset.site_pos_vs_ref_north ** 2), unit="meter"
        )
        dset.add_float("vpe", val=np.absolute(dset.site_pos_vs_ref_vertical), unit="meter")

    # Change sign of "gnss_satellite_clock" field
    if "delay.gnss_satellite_clock" in dset.fields:
        dset.delay.gnss_satellite_clock[:] = -dset.delay.gnss_satellite_clock

    # Change sign of "residual" field
    if "residual" in dset.fields:
        
        if pipeline == "gnss_vel":
            dset.residual[:] = -dset["doppler residual"]
        else:
            dset.residual[:] = -dset.residual
            
    if "APO" in file_.name:
        
        ref_pos_def = { "krss": [3348188.52082, 465041.62228, 5390743.08628],
                        "koug": [3855263.34560, -5049731.99430, 563040.40140] }
        ref_pos = position.Position(val=ref_pos_def[station], system="trs")
        enu = (dset.site_pos.trs.pos - ref_pos).enu
        dset.add_float("site_pos_vs_ref_east", val=enu.east, unit="meter")
        dset.add_float("site_pos_vs_ref_north", val=enu.north, unit="meter")
        dset.add_float("site_pos_vs_ref_up", val=enu.up, unit="meter")
        
        # Add HPE and VPE to dataset
        dset.add_float("hpe", val=np.sqrt(enu.east ** 2 + enu.north ** 2), unit="meter")
        dset.add_float("vpe", val=np.absolute(enu.up), unit="meter")

    # Write dataset
    dset.write_as(rundate=rundate, pipeline=pipeline, stage=stage, station=station)

    return dset


############################################################################################
#
#
#                                MAIN PROGRAM
#
#
############################################################################################
if __name__ == "__main__":
    """

    ADOP<wwwwd>_0000.csv - DOP/PL
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | Constellation |      | Constellation number (1:GPS, 2:GLONASS, 3:SBAS, 4:GALILEO, 5:BEIDOU, 6:QZSS)        |
    | SatInView     |      | Number of satellites in view                                                        |
    | UsedSat       |      | Number of used satellites                                                           |
    | HDOP          |      | HDOP                                                                                |
    | VDOP          |      | VDOP                                                                                |
    | PDOP          |      | PDOP                                                                                |
    | TDOP          |      | TDOP                                                                                |
    | GDOP          |      | GDOP                                                                                |
    | HPE           |      | HPE                                                                                 |
    | VPE           |      | VPE                                                                                 |

    AISD<wwwwd>_0000.csv - Ionospheric Delays
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | PRN           |      | Satellite number                                                                    |
    | Azimuth       |  deg | Azimuth (degrees clockwise from north direction)                                    |
    | Elevation     |  deg | Elevation (degrees from horizontal)                                                 |
    | KLBC          |    m | Ionospheric delay of Klobuchar model                                                |
    | NeQuick       |    m | Ionospheric delay of Nequick model                                                  |

    ALOS<wwwwd>_0000.csv - Satellite data
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | PRN           |      | Satellite number                                                                    |
    | Azimuth       |  deg | Azimuth (degrees clockwise from north direction)                                    |
    | Elevation     |  deg | Elevation (degrees from horizontal)                                                 |
    | LocalX        |      | Local vector (X)                                                                    |
    | LocalY        |      | Local vector (Y)                                                                    |
    | LocalZ        |      | Local vector (Z)                                                                    |
    | XCorrection   |    m | SV position correction (X)                                                          |
    | YCorrection   |    m | SV position correction (Y)                                                          |
    | ZCorrection   |    m | SV position correction (Z)                                                          |

    APO<wwwwd>_0000.csv - Positioning data
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | Latitude      |  deg | WGS84 latitude (degrees)                                                            |
    | Longitude     |  deg | WGS84 longitude (degrees)                                                           |
    | Height        |  deg | WGS84 height (degrees)                                                              |

    APRP<wwwwd>_0000.csv - Pseudo-range processing
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | PRN           |      | Satellite number                                                                    |
    | Signal        |  deg | Signal number                                                                       |
    | PseudoRange   |    m | Raw pseudo range                                                                    |
    | Phase         |  cyc | Raw phase measurement                                                               |
    | Doppler       |   Hz | Doppler measurement                                                                 |
    | SmoothedPseudoRange   |    m | Smoothed pseudo range                                                       |
    | SmoothingFilterWeight |    m | Weight of the smoothing filter                                              |
    | GroupDelay    |    m | Group delay                                                                         |
    | ISC           |    m | Inter-signal correction                                                             |
    | SignalStrength| dBHz | Signal strength                                                                     |
    | Frequency     |   Hz | Signal frequency                                                                    |

    ARES<wwwwd>_0000.csv - Tropospheric delays and residuals
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | PRN           |      | Satellite number                                                                    |
    | Signal        |  deg | Signal number                                                                       |
    | Azimuth       |  deg | Azimuth (degrees clockwise from north direction)                                    |
    | Elevation     |  deg | Elevation (degrees from horizontal)                                                 |
    | Residual      |    m | Residuals                                                                           |
    | Residual (Absolute Value) |   m | Residuals in absolute value                                              |
    | Doppler Residual |   m | Doppler residuals                                                                 |
    | Doppler Residual (Absolute Value)|   m | Doppler residuals in absolute value                               |
    | TropoDelay    |    m | Tropospheric delay                                                                  |

    ARPE<wwwwd>_0000.csv - Position errors vs reference
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | 3DvsRef       |    m | 3D position error vs reference                                                      |
    | 2DvsRef       |    m | Horizontal position error vs reference                                              |
    | XvsRef        |    m | X position error vs reference                                                       |
    | YvsRef        |    m | Y position error vs reference                                                       |
    | ZvsRef        |    m | Z position error vs reference                                                       |
    | VerticalvsRef |    m | Vertical position error vs reference                                                |
    | LatitudevsRef |  deg | Latitude position error vs reference                                                |
    | LongitudevsRef |  deg | Longitude position error vs reference                                              |
    | NorthvsRef    |    m | North position error vs reference                                                   |
    | EastvsRef     |    m | East position error vs reference                                                    |

    ASIV<wwwwd>_0000.csv - Satellites in view per constellation
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | Constellation |      | Constellation number (1:GPS, 2:GLONASS, 3:SBAS, 4:GALILEO, 5:BEIDOU, 6:QZSS)        |
    | SatInView     |      | Number of satellites in view per constellation                                      |
    | UsedSat       |      | Number of used satellites per constellation                                         |

    ASPO<wwwwd>_0000.csv - Satellite position/clock
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | PRN           |      | Satellite number                                                                    |
    | XPos          |    m | WGS84 X position                                                                    |
    | YPos          |    m | WGS84 Y position                                                                    |
    | ZPos          |    m | WGS84 Z position                                                                    |
    | Latitude      |  deg | WGS84 latitude (degrees)                                                            |
    | Longitude     |  deg | WGS84 longitude (degrees)                                                           |
    | Height        |    m | WGS84 height (degrees)                                                              |
    | Clock         |    m | Clock offset                                                                        |
    | IODC          |      | IODC                                                                                |
    | IODE          |      | IODE                                                                                |
    | EphemStatus   |      | Status of the satellite from ephemeris                                              |

    ATIM<wwwwd>_0000.csv - Time data
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | Clock Offset  |    s | Receiver clock offset                                                               |
    | Clock Drift   | ns/s | Receiver clock drift                                                                |
    | Constellation |      | Constellation number (1:GPS, 2:GLONASS, 3:SBAS, 4:GALILEO, 5:BEIDOU, 6:QZSS)        |

    AVEL<wwwwd>_0000.csv - Velocity data
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | 3DSpeed       |  m/s | Instantaneous 3D speed                                                              |
    | HSpeed        |  m/s | Instantaneous horizontal speed                                                      |
    | XSpeed        |  m/s | Instantaneous speed along WGS84 X axis                                              |
    | YSpeed        |  m/s | Instantaneous speed along WGS84 Y axis                                              |
    | ZSpeed        |  m/s | Instantaneous speed along WGS84 Z axis                                              |
    | NorthSpeed    |  m/s | Instantaneous north speed                                                           |
    | EastSpeed     |  m/s | Instantaneous east speed                                                            |
    | VerticalSpeed |  m/s | Instantaneous vertical speed                                                        |
    | Heading       |  deg | Instantaneous heading speed                                                         |

    RNAV<wwwwd>_0000.csv - Navigation Parameters
    | Field         | Unit | Description                                                                         |
    |---------------|------|-------------------------------------------------------------------------------------|
    | PRN           |      | Satellite number                                                                    |
    | TOC           |      | Time of clock                                                                       |
    | TOE           |      | Time of ephemeris                                                                   |
    | IODC          |      | IODC                                                                                |
    | IODE          |      | IODE                                                                                |
    | SV Accuracy   |      | SV Accuracy                                                                         |
    | SV Health     |      | SV Health                                                                           |
    | TGD           |      | Group Delay                                                                         |
    """
    dset_diff_merged = None
    log.init(log_level="debug")
    #year = 2019
    #month = 7
    #day = 1
    year = 2020
    month = 1
    day = 4
    # wwwwd = "20473"
    # wwwwd = "20556"
    #wwwwd = "20601"  # 2019-07-01
    #wwwwd = "20863" # 2020-01-01
    wwwwd = "20866" # 2020-01-04
    rundate = date(year, month, day)
    #station = "nabf"
    #station = "krss"
    station = "koug"
    stage = "write"
    pipeline = "gnss_vel"  # gnss, gnss_vel
    system = "E"  # G, E

    # directory = pathlib.Path("./data/HOFS_093_KLOB_L1_GAL/20190403")
    # directory = pathlib.Path("./data/vegs_152_KLOB_L1_GAL/20190601")
    #directory = pathlib.Path(f"./data/KRSS_001_SPV_NEQ_L1_GAL/20200101")
    #directory = pathlib.Path(f"./data/NABD_001_SPV_NEQ_L1_GAL/20200101")
    directory = pathlib.Path(f"./data/KOUG_004_SPV_NEQ_L1_GAL/20200104")

    #directory = pathlib.Path(f"./data/{station}_182_KLOB_L1_GPS/20190701")    
    #directory = pathlib.Path(f"./data/{station}_182_L1L2_GPS/20190701") 
    #directory = pathlib.Path(f"./data/{station}_182_NEQ_L1_GAL/20190701")
    #directory = pathlib.Path(f"./data/{station}_182_KLOB_L1_GAL/20190701")
    #directory = pathlib.Path(f"./data/{station}_182_E1E5A_GAL/20190701")
    #directory = pathlib.Path(f"./data/{station}_182_E1E5B_GAL/20190701")

    file_overview = [
        # f"ADOP{wwwwd}_0000.csv",  # DOP/PL
        # f"AIPD{wwwwd}_0000.csv",  # Ionospheric Pierce Point data
        # f"AISD{wwwwd}_0000.csv",  # Ionospheric Delays
        # f"ALOS{wwwwd}_0000.csv",  # Satellite data
        # f"AMPE{wwwwd}_0000.csv",  # Position errors vs mean
        #f"APO{wwwwd}_0000.csv",  # Positioning data
        #f"APRP{wwwwd}_0000.csv",  # Pseudo-range processing data
        f"ARES{wwwwd}_0000.csv",  # Tropospheric delays and residuals
        #f"ARPE{wwwwd}_0000.csv",  # Position errors vs reference
        # f"ARVE{wwwwd}_0000.csv",  # Velocity data
        # f"ASAS{wwwwd}_0000.csv",  # Service Availability Statistics
        # f"ASIV{wwwwd}_0000.csv",  # Satellites in view per constellation
        #f"ASPO{wwwwd}_0000.csv",  # Position/clock
        # f"ATIM{wwwwd}_0000.csv",  # Time data
        # f"AURE{wwwwd}_0000.csv",  # Tropospheric/noise/multipath data
        #f"AVEL{wwwwd}_0000.csv",  # Velocity errors vs reference
        #f"RNAV{wwwwd}_0000.csv",  # Navigation Parameters
    ]

    # Initialize configuration
    config.init(rundate, pipeline, station=station)

    # Where dataset
    dset_where = _generate_where_dset(rundate, pipeline, station, stage, system)

    files = [directory / file_ for file_ in file_overview]
    for file_ in files:

        # Read Spring output files
        dset_spring = _read_spring_files(file_, rundate, pipeline, station)

        ## Generate difference
        #index_by = "time, satellite" if "satellite" in dset_where.fields else "time"
        #dset_diff = dset_where.difference(dset_spring, index_by=index_by)

        ## Merge differenced datasets
        #if dset_diff_merged is None:
        #    dset_diff_merged = dset_diff
        #else:
        #    _merge_diff_dset(dset_diff_merged, dset_diff)

    # Save differences
    #dset_diff_merged.analysis.update(config.date_vars(rundate))
    #dset_diff_merged.write_as(pipeline=pipeline, stage="diff_where_spring", station=station)
