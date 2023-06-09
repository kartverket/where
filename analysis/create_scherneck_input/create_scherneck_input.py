import sys

from where import apriori
from where.lib import config
from where.data.time import Time

"""
Prints site_id and X,Y,Z values to the console. Format compatible with the form at http://holt.oso.chalmers.se/loading/.

Note! The form only accept 100 stations at the time. Example of usage:

python3 create_scherneck_input.py vlbi > out
python3 create_scherneck_input.py slr >> out
split -l100 out scherneck_input

This gives files with maximum 100 stations each:
scherneck_inputaa  scherneck_inputab  scherneck_inputac
"""


def main():
    try:
        tech = sys.argv[1]
    except IndexError:
        print("Missing argument <tech>")
        sys.exit(1)

    if tech == "vlbi":
        reference_frames = "itrf, vtrf"
    else:
        reference_frames = "itrf, custom"

    config.set_analysis(rundate=None, pipeline=tech)
    config.files.profiles = [tech]
    config.read_pipeline(tech)

    trf = apriori.get("trf", time=Time.now(), reference_frames=reference_frames)

    out_fmt = "{:<24} {:>16.3f}{:>16.3f}{:>16.3f}"

    for site in trf:
        print(out_fmt.format(site.key, site.pos.trs.x, site.pos.trs.y, site.pos.trs.z))


if __name__ == "__main__":
    main()
