# Standard library imports
import argparse
from datetime import datetime
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
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("station", help="Station name or two letter code", type=str)
parser.add_argument("--start_year", help="Start year", type=int)
parser.add_argument("--end_year", help="End year", type=int)

# Setup config to make apriori modules work
pipeline = "vlbi"
config.set_analysis(rundate=None, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

args = parser.parse_args()
start = args.start_year if args.start_year else 1979 
end = args.end_yaer if args.end_year else datetime.now().year
station = args.station

# Create output directories
csvdir = "csv"
os.makedirs(csvdir, exist_ok=True)
csvfile = f"{csvdir}/qcodes_{station}_{start}_{end}.csv"

QCODES = [ '0',  '1',  '2',  '3',  '4',  '5',  '6',  '7',  '8', '9',  'G',  'H',  'N',  '-']
header = ["VgosDB"] + [f"Total_{code}" for code in QCODES]
header = ",".join(header) + "\n"

with open(csvfile, 'w') as fid:
    fid.write(header)
    for year in range(start, end + 1):
        master = apriori.get("vlbi_master_schedule", rundate=datetime(year,1, 1))
        for session in master:
            rundate = datetime.strptime(master[session]["date"],"%Y%m%d")
            file_vars = config.create_file_vars(rundate, pipeline, session_code=session)
            config.init(rundate, pipeline, session_code=session)
            pipeline_file_vars = plugins.call(
                package_name="where.pipelines", plugin_name=pipeline, part="file_vars", file_vars=file_vars
            )
            file_vars.update(pipeline_file_vars)
            parser = parsers.parse_key(f"vlbi_obs_vgosdb", file_vars)
            report = parser.data["correlation_report"]
            
            sta_letter = report["stations"][station]
            sta_idx = np.char.find(report["qcodes"]["bl"], sta_letter) > 0
            sum_qcodes = [str(np.sum(report["qcodes"][code][sta_idx])) for code in QCODES]
            sum_qcodes.insert(0, file_vars["input_data_name"])
            line = ','.join(sum_qcodes) + "\n"
            fid.write(line) 
