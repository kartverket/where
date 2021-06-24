"""Write SISRE analysis results

Description:
------------
TODO: Has IGS monitoring group defined a SISRE format?


"""
# Standard library imports
from collections import namedtuple

# External library imports
import numpy as np

# Midgard imports
from midgard.dev import console
from midgard.dev import plugins
from midgard.writers._writers import get_field, get_header

# Where imports
import where
from where.lib import config
from where.lib import util
from where import pipelines
from where.writers import sisre_output_buffer

WriterField = namedtuple(
    "WriterField", ["name", "field", "attrs", "dtype", "format", "width", "header", "unit", "description"]
)
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        name  (str):             Unique field name
        field (str):             Dataset field name
        attrs (Tuple[str]):      Field attributes
        dtype (Numpy dtype):     Type of field
        format (str):            Format string
        width (int):             Width of header information
        header (str):            Header information
        unit (str):              Unit of field
        description (str):       Description of field
    """

FIELDS = (
    WriterField(
        "date",
        "date",
        (),
        object,
        "%21s",
        19,
        "DATE",
        "YYYY/MM/DD hh:mm:ss",
        "Date in format year, month, day and hour, minute and second",
    ),
    WriterField("mjd", "time", ("gps", "mjd"), float, "%14.6f", 14, "MJD", "", "Modified Julian Day"),
    WriterField("gpsweek", "time", ("gps", "gps_ws", "week"), int, "%5d", 5, "WEEK", "", "GPS week"),
    WriterField(
        "gpssec", "time", ("gps", "gps_ws", "seconds"), float, "%11.3f", 11, "GPSSEC", "second", "GPS seconds"
    ),
    WriterField("satellite", "satellite", (), object, "%5s", 5, "SAT", " ", "Satellite number"),
    WriterField(
        "used_iode",
        "used_iode",
        (),
        float,
        "%6d",
        6,
        "IODE",
        "",
        f"Used ephemeris issue of data indicates changes to the broadcast ephemeris:\n"
        f"""
{'': >38}- GPS:     Ephemeris issue of data (IODE), which is set equal to IODC
{'': >38}- Galileo: Issue of Data of the NAV batch (IODnav)
{'': >38}- QZSS:    Ephemeris issue of data (IODE)
{'': >38}- BeiDou:  Age of Data Ephemeris (AODE)
""",
    ),
    WriterField(
        "trans_time_gpsweek",
        "trans_time_gpsweek",
        (),
        object,
        "%15s",
        15,
        "TRANS_TIME",
        "wwwwd:ssssss",
        "Transmission time (receiver reception time)",
    ),
    WriterField("toe_gpsweek", "toe_gpsweek", (), object, "%15s", 15, "TOE", "wwwwd:ssssss", "Time of ephemeris"),
    WriterField(
        "diff_trans_toe",
        "diff_trans_toe",
        (),
        float,
        "%8d",
        8,
        "TM-TOE",
        "second",
        "Difference between transmission time and time of ephemeris",
    ),
    WriterField(
        "age_of_ephemeris",
        "age_of_ephemeris",
        (),
        float,
        "%8d",
        8,
        "T-TOE",
        "second",
        "Age of ephemeris, which is the difference between the observation time and the time of ephemeris " "(ToE)",
    ),
    # WriterField("diff_time_trans", "diff_time_trans", (), float, "%8d", 8, "T-TM", "", ""),
    WriterField(
        "clk_diff",
        "clk_diff",
        (),
        float,
        "%16.4f",
        16,
        "ΔCLOCK",
        "meter",
        "Satellite clock correction difference related to center of mass of satellite and corrected for "
        "satellite biases",
    ),
    WriterField(
        "clk_diff_with_dt_mean",
        "clk_diff_with_dt_mean",
        (),
        float,
        "%16.4f",
        16,
        "ΔCLOCK_MEAN",
        "meter",
        "Satellite clock correction difference related to center of mass of satellite and corrected for "
        "satellite biases and averaged clock offset in each epoch",
    ),
    WriterField(
        "dalong_track",
        "dalong_track",
        (),
        float,
        "%16.4f",
        16,
        "ΔALONG_TRACK",
        "meter",
        "Satellite coordinate difference between broadcast and precise ephemeris in along-track direction",
    ),
    WriterField(
        "dcross_track",
        "dcross_track",
        (),
        float,
        "%16.4f",
        16,
        "ΔCROSS_TRACK",
        "meter",
        "Satellite coordinate difference between broadcast and precise ephemeris in cross-track direction",
    ),
    WriterField(
        "dradial",
        "dradial",
        (),
        float,
        "%16.4f",
        16,
        "ΔRADIAL",
        "meter",
        "Satellite coordinate difference between broadcast and precise ephemeris in radial direction",
    ),
    WriterField("orb_diff_3d", "orb_diff_3d", (), float, "%16.4f", 16, "ORB_DIFF_3D", "meter", "3D orbit difference"),
    WriterField("sisre_orb", "sisre_orb", (), float, "%16.4f", 16, "SISRE_ORB", "meter", "Orbit-only SISRE"),
    WriterField("sisre", "sisre", (), float, "%16.4f", 16, "SISRE", "meter", "Global averaged SISRE"),
)


@plugins.register
def sisre_writer(dset: "Dataset") -> None:
    """Write SISRE analysis results

    Args:
        dset:   A dataset containing the data.
    """
    # Add additional fields used by the writer
    dset.add_text("date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])
    dset.add_text(
        "time_gpsweek", 
        val=[f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.time],
        write_level="detail",
    )
    dset.add_text(
        "trans_time_gpsweek",
        val=[
            f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.used_transmission_time
        ],
        write_level="detail",
    )
    dset.add_text(
        "toe_gpsweek",
        val=[f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.used_toe],
        write_level="detail",
    )
    # dset.add_float("diff_time_trans", val=(dset.time.mjd - dset.used_transmission_time.mjd) * Unit.day2second, Unit="second")
    dset.add_float(
        "dalong_track", 
        val=dset.orb_diff.acr.along, 
        unit=dset.unit("orb_diff.acr.along"), 
        write_level="detail",
    )
    dset.add_float(
        "dcross_track", 
        val=dset.orb_diff.acr.cross, 
        unit=dset.unit("orb_diff.acr.cross"), 
        write_level="detail",
    )
    dset.add_float(
        "dradial", 
        val=dset.orb_diff.acr.radial, 
        unit=dset.unit("orb_diff.acr.radial"),
        write_level="detail",
    )

    ## Add 'detail' fields used by the writer
    # write_level = config.tech.get("write_level", default="operational").as_enum("write_level")
    # if write_level <= enums.get_value("write_level", "detail"):
    #    FIELDS += (
    #        WriterField("clk_brdc_com", "clk_brdc_com", (), float, "%16.4f", 16, "CLK_BRDC", ""),
    #        WriterField("clk_precise_com", "clk_precise_com", (), float, "%16.4f", 16, "CLK_PRECISE", ""),
    #        WriterField("bias_brdc", "bias_brdc", (), float, "%10.4f", 10, "B_BRDC", ""),
    #        WriterField("bias_precise", "bias_precise", (), float, "%10.4f", 10, "B_PREC", ""),
    #        WriterField("dt_mean", "dt_mean", (), float, "%10.4f", 10, "dt_MEAN", ""),
    #    )
    #
    #    dset.add_float("dt_mean", val=dset.clk_diff - dset.clk_diff_with_dt_mean, unit="meter", write_level="detail")

    # List epochs ordered by satellites
    idx = np.concatenate([np.where(dset.filter(satellite=s))[0] for s in dset.unique("satellite")])

    # Put together fields in an array as specified by the fields-tuple
    output_list = list(zip(*(get_field(dset, f.field, f.attrs, f.unit) for f in FIELDS)))
    output_array = np.array(output_list, dtype=[(f.name, f.dtype) for f in FIELDS])[idx]

    # Write to disk
    # NOTE: np.savetxt is used instead of having a loop over all observation epochs, because the performance is better.
    file_path = config.files.path(f"output_sisre_{dset.vars['label']}", file_vars={**dset.vars, **dset.analysis})
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in FIELDS),
        header=_get_header(dset),
        delimiter="",
        encoding="utf8",
    )

    # Append SISRE output path to SISRE output buffer file
    if config.tech.sisre_writer.write_buffer_file.bool:
        sisre_output_buffer.sisre_output_buffer(dset)


def _get_header(dset: "Dataset") -> str:
    """Get header

    Args:
        dset:   A dataset containing the data.

    Returns:
        Header lines
    """
    # SISRE configuration
    add_description = "\nSISRE ANALYSIS CONFIGURATION\n\n"
    add_description += str(config.tech.as_str(key_width=25, width=70, only_used=True)) + "\n\n\n"
    add_description += _get_paths()

    # Information about used biases and phase center offsets (PCOs)
    add_description += f"{'SAT':^4s}{'BIAS_BRDC':>10s}{'BIAS_PREC':>10s}{'PCO_BRDC':^26s}{'PCO_PREC':^26s}\n"
    add_description += f"{'':^4s}{'[m]':^10s}{'[m]':^10s}{'[m]':^26s}{'[m]':^26s}\n"
    for sat in dset.unique("satellite"):
        add_description += "{:>3s}{:>10.4f}{:>10.4f}{:>10.4f}{:>8.4f}{:>8.4f}{:>10.4f}{:>8.4f}{:>8.4f}\n".format(
            sat,
            dset.meta["bias_brdc"][sat],
            dset.meta["bias_precise"][sat],
            dset.meta["pco_sat_brdc"][sat][0],
            dset.meta["pco_sat_brdc"][sat][1],
            dset.meta["pco_sat_brdc"][sat][2],
            dset.meta["pco_sat_precise"][sat][0],
            dset.meta["pco_sat_precise"][sat][1],
            dset.meta["pco_sat_precise"][sat][2],
        )

    header = get_header(
        FIELDS,
        pgm_version=f"where {where.__version__}",
        run_by=util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else "",
        summary="SISRE analysis results",
        add_description=add_description + "\n\n",
    )

    return header


def _get_paths() -> str:
    """Get file paths of used files

    Returns:
        Header file path lines
    """
    lines = "[file_paths]\n"
    key_width = 25
    fill_args = dict(width=120, hanging=key_width + 3, break_long_words=False, break_on_hyphens=False)
    path_def = {
        "Broadcast orbit": "gnss_rinex_nav_.",
        "Precise orbit": "gnss_orbit_sp3",
        "Bias": "gnss_sinex_bias",
        "Precise clock": "gnss_rinex_clk",
    }

    for name, file_key in sorted(path_def.items()):
        if file_key == "gnss_rinex_clk":
            if config.tech.clock_product.str != "clk":
                continue
        path = ", ".join(str(p) for p in sorted(pipelines.paths(file_key)))
        lines += console.fill(f"{name:<{key_width}} = {path}", **fill_args)
        lines += "\n"

    lines += "\n\n"
    return lines
