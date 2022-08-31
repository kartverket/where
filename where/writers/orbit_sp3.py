"""Write SP3-d orbit file

Description:
------------
Write data in SP3-d file format (see :cite:`hilla2016`). The reference frame is defined in the header of the SP3 file,
which is for IGS products the current IGS terrestrial reference frame (e.g. IGb08). The time system is for IGS products
the GPS time scale. 

Example:
--------

from where import writers
writers.write(dset)


"""
# Standard library imports
from math import floor
from typing import List, Set, Union

# Midgard imports
from midgard.data.time import Time
from midgard.dev import plugins
from midgard.math.constant import constant
from midgard.math.unit import Unit

# Where imports
from where.lib import config

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])

@plugins.register
def orbit_sp3(dset: "Dataset"):
    """Write SP3-d orbit file

    Args:
        dset:  Dataset, a dataset containing the data.
    """

    time = Time(
            val=dset.time.gps.datetime[0],
            scale="gps", 
            fmt="datetime",
    )

    with config.files.open("output_orbit_sp3", file_vars={**dset.vars, **dset.analysis}, mode="wt") as fid:

        #
        # Write SP3 header
        #

        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # #dP2016  3  1  0  0  0.00000000      96 ORBIT IGb08 HLM  IGS
        fid.write(
            "{:>2s}{:1s}{:4.0f} {:2.0f} {:2.0f} {:2.0f} {:2.0f} {:11.8f} {:7d} {:5s} {:5s} {:3s} {:>4s}\n"
            "".format(
                "#d", # Version
                "P", # Pos or Vel flag
                time.year, # Year start
                time.month, # Month start
                time.day, # Day start
                time.hour, # Hour start
                time.minute, # Minute start
                time.second, # Second start
                len(dset.unique("time")),
                "ORBIT", # Data used
                "GTRF", # Coordinate system
                "BCT", # Orbit type
                "NMA", # Agency
            )
        )

        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # ## 1886 172800.00000000   900.00000000 57448 0.0000000000000
        fid.write(
            "## {:4.0f} {:15.8f} {:14.8f} {:5.0f} {:15.13f}\n"
            "".format(
                time.gps_ws.week, # GPS week
                time.gps_ws.seconds, # Seconds of week
                dset.meta["sampling_rate"], # Epoch interval
                time.mjd_int, # Modified Julian day - integer part
                time.mjd_frac, # Modified Julian day - fractional part
            )
        )

        # ----+----1----+----2----+----3----+----4----+----5----+----6        
        # +  116   G01G02G03G04G05G06G07G08G09G10G11G12G13G14G15G16G17
        # +        G18G19G20G21G22G23G24G25G26G27G28G29G30G31G32R01R02
        # +        R03R04R05R07R08R09R12R13R14R15R17R18R19R20R21R22R24
        # +        E01E02E03E04E05E07E08E10E11E12E13E14E15E18E19E21E24
        # +        E25E26E27E30E31E33E34E36C06C07C08C09C10C11C12C13C14
        # +        C16C19C20C21C22C23C24C25C26C27C28C29C30C32C33C34C35
        # +        C36C37C38C39C40C41C42C43C44C45C46J01J02J03  0  0  0
        satellites = sorted(set(dset.satellite))
        lines = _newline_satellite(satellites, 17)
        for idx, line in enumerate(lines):
            fid.write(
                "+  {:>3s}   {:51s}\n"
                "".format(
                    str(len(satellites)) if idx==0 else "   ",  # Number of satellites or empty string
                    line, 
                )
            )

        # At least 5 satellite lines has to be written
        for idx in range(0, 5 - len(lines)):
            fid.write("+        {:51s}\n".format("  0" * 17))

        # ----+----1----+----2----+----3----+----4----+----5----+----6         
        # ++         5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
        # ++         5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
        # ++         5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
        # ++         5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
        # ++         5  5  5  5  5  5  5  5  5  5  5  7  5  5  5  5  5
        # ++         5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5  5
        # ++         5  5  5  5  5  5  5  5  5  5  5  5  5  5  0  0  0
        accuracy = config.tech[_SECTION].accuracy.str
        accuracy = accuracy.rjust(3) if accuracy else "  0" # Default: Unknown. If accuracy is known, then it 

        lines = _newline_satellite([accuracy] * len(satellites), 17)
        for idx, line in enumerate(lines):
            fid.write(
                "++       {:51s}\n"
                "".format(
                    line,  
                )
            )

        # At least 5 accuracy lines has to be written
        for idx in range(0, 5 - len(lines)):
            fid.write("++       {:51s}\n".format("  0" * 17))
                
        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # %c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
        systems = dset.unique("system")
        fid.write(
            "%c {:2} cc {:3s} ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n"
            "".format(
                systems[0] if len(systems) == 1 else "M", # File type
                "GPS", # Time system 
            )
        )
        
        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
        fid.write("%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc\n")

        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # %f  1.2500000  1.025000000  0.00000000000  0.000000000000000
        fid.write(
            "%f {:10.7f} {:12.9f} {:14.11f} {:18.15f}\n"
            "".format(
                1.25, # Base for Pos/Vel (mm or 10^-4 mm/sec)
                1.025, # Base for Clk/Rate (psec or 10^-4 psec/sec)
                0.0, # 14-column float
                0.0, # 18-column float
            )
        )
        
        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # %f 0.0000000 0.000000000 0.00000000000 0.000000000000000
        fid.write(
            "%f {:10.7f} {:12.9f} {:14.11f} {:18.15f}\n"
            "".format(
                0.0, # 10-column float
                0.0, # 12-column float
                0.0, # 14-column float
                0.0, # 18-column float
            )
        )
        
        # ----+----1----+----2----+----3----+----4----+----5----+----6
        # %i    0    0    0    0      0      0      0      0         0
        # %i    0    0    0    0      0      0      0      0         0
        fid.write(
            "%i {:4d} {:4d} {:4d} {:4d} {:6d} {:6d} {:6d} {:6d} {:9d}\n"
            "".format(
                0, # 4-column int
                0, # 4-column int
                0, # 4-column int
                0, # 4-column int
                0, # 6-column int
                0, # 6-column int
                0, # 6-column int
                0, # 6-column int
                0, # 9-column int
            )
        )
        fid.write(
            "%i {:4d} {:4d} {:4d} {:4d} {:6d} {:6d} {:6d} {:6d} {:9d}\n"
            "".format(
                0, # 4-column int
                0, # 4-column int
                0, # 4-column int
                0, # 4-column int
                0, # 6-column int
                0, # 6-column int
                0, # 6-column int
                0, # 6-column int
                0, # 9-column int
            )
        )
        
        #
        # Write SP3 data block
        #
        # *  2016  3  1  0  0  0.00000000
        # PG01  10138.887745 -20456.557725 -13455.830128     13.095853  7  6  4 137
        # PG02 -21691.921884  13338.131173  -6326.904893    599.417359 10  7  8 114
        # PG03   1061.483783 -15622.426751 -21452.532447    -30.693182  7  6 10 125
        # PG04  25398.213954   6966.881030   4188.487313 999999.999999
        # PG05 -20431.536257   5143.439387  16220.688245   -141.634044  9  8  7 101
        # PG06 -19896.425641    760.039334 -17564.616810    153.322562 12 10 13 104
        # PG07  -5499.233721 -14614.719047  21564.675152    467.390140  7 10  7 123
        # PG08   7359.468852 -20268.667422  15550.334189    -25.778533  8  6  8 113
        for time in sorted(dset.unique("time")):
            idx = dset.filter(time=time)
            
            # Epoch header
            # ----+----1----+----2----+----3----+----4----+----5----+----6
            # *  2016  3  1  0  0  0.00000000
            fid.write(
                "*  {:4d} {:2d} {:2d} {:2d} {:2d} {:11.8f}\n"
                "".format(
                    dset.time.gps.year[idx][0], # Year start
                    dset.time.gps.month[idx][0], # Month start
                    dset.time.gps.day[idx][0], # Day start
                    dset.time.gps.hour[idx][0], # Hour start
                    dset.time.gps.minute[idx][0], # Minute start
                    dset.time.gps.second[idx][0], # Second start
                )
            )
            
            for sat in sorted(dset.unique("satellite")):
                
                if sat not in set(dset.satellite[idx]):
                    continue

                idx_sat = dset.satellite[idx] == sat
                     
                # Position and clock record
                # ----+----1----+----2----+----3----+----4----+----5----+----6
                # PG01  10138.887745 -20456.557725 -13455.830128     13.095853 
                # PG02 -21691.921884  13338.131173  -6326.904893    599.417359 
                fid.write(
                    "P{:3s}{:14.6f}{:14.6f}{:14.6f}{:14.6f}\n"
                    "".format(
                        dset.satellite[idx][idx_sat][0], # Vehicle identifier
                        dset.sat_posvel.pos.trs.x[idx][idx_sat][0] * Unit.meter2kilometer, # X-coordinate (km)
                        dset.sat_posvel.pos.trs.y[idx][idx_sat][0] * Unit.meter2kilometer, # Y-coordinate (km)
                        dset.sat_posvel.pos.trs.z[idx][idx_sat][0] * Unit.meter2kilometer, # Z-coordinate (km)
                        (-dset.delay.gnss_satellite_clock[idx][idx_sat][0] / constant.c) * Unit.second2microsecond, # Clock (microsec)
                    )
                )


        # End of file
        # ----+----1----+----2----+----3----+----4----+----5----+----6
        fid.write("EOF\n")
        
def _newline_satellite(words: Union[List[str], Set[str]], num_words: int) -> List[str]:
    """Generate a string, whereby a newline is set after a defined number of words

    Args:
        words:     List with words
        num_words: Number of words after which a newline should be set

    Returns:
        String with newlines
    """
    words = list(words)
    newlines = list()
    for idx in range(0, len(words), num_words):
        line = "".join(words[idx : idx + num_words])
        
        # Fill line length with 0 after all satellites identifiers are listed
        if not len(line) == num_words*len(words[0]):  
            diff = num_words*len(words[0]) - len(line)
            line = line + "  0" * floor(diff/3)
            
        newlines.append(line)
    return newlines       
        
        
