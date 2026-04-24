""" Reads vlbi correlator reports and creates csv files for each station


QCODES definitions
* 0        no fringe detected

* 1-9      fringe detected, higher value means better quality

* B        fourfit interpolation error

* D        no data in one or more frequency channels

* E        fringe found at edge of SBD, MBD, or rate window

* F        fork problem in processing

* G        channel amplitude diverges too far from mean amplitude

* H        low phase-cal amplitude in one or more channels

* N        correlation or fringing failed

* -        correlation not attempted
"""

# Standard library imports
import argparse
from datetime import datetime
import glob
import pathlib
import os

# Third library imports
import numpy as np

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where import parsers
from where.lib import config

# Setup input argument parser for script 
parser = argparse.ArgumentParser(epilog="Example: python vlbi_qcodes.py Nn --start_year=2024 --end_year=2025",
                                 description="Reads vlbi correlator reports from vgosdb and accumulates the qcodes for each specified station. Results are written to a csv file.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("station", help="Station names or two letter codes", nargs="+", type=str)
parser.add_argument("--start_year", help="Read reports from this year and onwards", type=int, default=1979)
parser.add_argument("--end_year", help="Read reports up until and including this year", type=int, default=datetime.now().year)
parser.add_argument("--data_dir", help="Read vgosdb in specific folder. When not set the sessions from the master file will be searched for", type=str, default=None)

# Setup config to make apriori modules work
pipeline = "vlbi"
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

args = parser.parse_args()
start = args.start_year
end = args.end_year
stations = args.station
data_dir = pathlib.Path(args.data_dir) if args.data_dir else None

# Create output directories
csvdir = "csv"
os.makedirs(csvdir, exist_ok=True)

QCODES = [ '0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8', '9',  'B', 'D', 'E', 'F', 'G',  'H',  'N',  '-']
header = ["VgosDB"] + [f"Total_{code}" for code in QCODES]
header = ",".join(header) + "\n"


def process_file(vgosdb, file_path, stations, fids):
    """ Reads a correlation report and writes the output to a line in a csv file.
    """
    #print(f"Parsing {vgosdb}")
    parser = parsers.parse_file(f"vlbi_correlation_report", file_path=file_path)
    
    if not file_path.exists():
        print(f"No correlation report found. Skipping {vgosdb}")
        return
    report = parser.data
    if not report:
        print(f"No correlation report with supported format. Skipping {vgosdb}")
        return
    for s in stations:
        try:
            sta_letter = report["stations"][s]
        except KeyError:
            # station not in session
            return
        sta_idx = np.char.find(report["qcodes"]["bl"], sta_letter) >= 0
        sum_qcodes = [str(np.sum(report["qcodes"][code][sta_idx])) if code in report["qcodes"].dtype.names else '0' for code in QCODES]
        sum_qcodes.insert(0, vgosdb)
        line = ','.join(sum_qcodes) + "\n"
        fids[s].write(line)

fids = {}
for s in stations:
    csvfile = f"{csvdir}/qcodes_{s}_{start}_{end}.csv"
    fids[s] = open(csvfile, 'w')
    fids[s].write(header)

if data_dir:
    for file in glob.glob("**/*V000_kMk4.hist", root_dir=data_dir, recursive=True):
        file = pathlib.Path(file)
        vgosdb = file.stem.split("_")[0]
        file_path = data_dir / file
        process_file(vgosdb, file_path, stations, fids)
else:
    for year in range(start, end + 1):
        master = apriori.get("vlbi_master_schedule", rundate=datetime(year,1, 1))
        for session in master:
            rundate = datetime.strptime(master[session]["date"],"%Y%m%d")
            file_vars = config.create_file_vars(rundate, pipeline, session_code=session)
            vgosdb = f"{master[session]['date']}-{session}"
            file_path = config.files.path("vlbi_corr_report_vgosdb", file_vars)
            process_file(vgosdb, file_path, stations, fids)

for fid in fids.values():
    fid.close()
