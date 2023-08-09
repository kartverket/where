import argparse
import subprocess
import glob
import sys
import datetime

from midgard.config import config as mg_config

from where.lib import config
from where import parsers

def run_session(session, setup, write_to_library):
    cmd = ["where"]
    cmd += session
    setup_name = setup.split(".")[0]
    cfg = mg_config.Configuration.read_from_file(setup_name, f"{setup}")
    for section, values in cfg.as_dict().items():
        for k, v in values.items():
            cmd.append(f"--{section}:{k}={v}")
    cmd = cmd + ["-v", "-N", f"--id={setup_name}", f"--write_to_library={write_to_library}"]
    subprocess.run(cmd)
    #print(cmd)

def create_session_list(stations, years, all_stations=True):
    session_list = []
    func = all if all_stations else any
    for year in years:
        yy = f"{year}"[-2:]
        data = parsers.parse_key("vlbi_master_file", file_vars=dict(yy=yy, yyyy=year)).as_dict()
        for key, session_dict in data.items():
            # ignore stations after the minus sign since these did not participate in the session after all
            stations_in_session = session_dict["stations"].split("-")[0]
            if func(station in stations_in_session for station in stations):
                session_date = datetime.datetime.strptime(f'{year} {session_dict["doy"]}', '%Y %j')
                if session_date < datetime.datetime.today():
                    session = [f"{year}", f"{session_date.month}", f"{session_date.day}", f"--session_code={key}"]
                    session_list.append(session)
    return session_list 

def main():

    parser = argparse.ArgumentParser(epilog="Example: python run_sessions.py --start 2020-01-01 --end 2021-02-15 --stations Ns Ny --all_stations --years 2020 2021 --setups nyale13s0 wettzell0")
    parser.add_argument("--start", help="Start date to look for sessions in master files. Format: YYYY-mm-dd", type=datetime.date.fromisoformat, default=datetime.date.min)
    parser.add_argument("--end", help="End date to look for sessions in master files. Format:YYYY-mm-dd", type=datetime.date.fromisoformat, default=datetime.date.max)
    parser.add_argument("--all_stations", help="Enable this flag if all stations in the list must be in the session. Otherwise only one stations in the list will be needed to process the session", action="store_true")
    parser.add_argument("--stations", help="Process sessions containing stations in this list.", nargs='+', default=["Ns", "Ny"])
    parser.add_argument("--years", help="Search for sessions in the master files from the years in this list.", default=[2020, 2021], nargs='+', type=int)
    parser.add_argument("--setups", help="Name of configurations to be run and dataset id to store it under. If omitted all existing configurations will be run", nargs='+')
    parser.add_argument("--write_to_library", help="Enable writing to config library", action="store_true")
    args = parser.parse_args()
    
    config.read_pipeline("vlbi")
    
    session_list = create_session_list(stations=args.stations, years=args.years, all_stations=args.all_stations)
    print(f"Found the following sessions: {session_list}")
    if args.setups is None:
        setups = glob.glob("*.conf")
    else:
        setups = [f"{setup}.conf" for setup in args.setups]
    print(f"Found the following {setups}")
    
    for session in session_list:
        date = datetime.date(int(session[0]), int(session[1]), int(session[2]))
        if date < args.start or date > args.end:
            continue

        for setup in setups:
            run_session(session, setup, args.write_to_library)

if __name__ == "__main__":
    main()
