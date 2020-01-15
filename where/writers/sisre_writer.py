"""Write SISRE analysis results

Description:
------------
TODO: Has IGS monitoring group defined a SISRE format?


"""
# Standard library imports
from collections import namedtuple
from datetime import datetime
from typing import Tuple

# External library imports
import numpy as np

# Midgard imports
import midgard
from midgard.dev import console
from midgard.dev import plugins

# Where imports
import where
from where.lib import config
from where.lib import util
from where import pipelines
from where.writers import sisre_output_buffer

WriterField = namedtuple("WriterField", ["field", "attrs", "dtype", "format", "width", "header", "unit"])
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        field (str):             Dataset field name
        attrs (Tuple[str]):      Field attributes
        dtype (Numpy dtype):     Type of field
        format (str):            Format string
        width (int):             Width of header information
        header (str):            Header information
        unit (str):              Unit of field
    """


@plugins.register
def sisre_writer(dset: "Dataset") -> None:
    """Write SISRE analysis results

    Args:
        dset:   A dataset containing the data.
    """
    write_level = config.tech.get("write_level", default="operational").as_enum("write_level")

    fields = (
        WriterField("time_date", (), object, "%21s", 19, "EPOCH", "YYYY/MM/DD hh:mm:ss"),
        WriterField("time", ("gps", "mjd"), float, "%14.6f", 14, "", "mjd"),
        WriterField("time_gpsweek", (), object, "%15s", 15, "", "wwwwd:ssssss"),
        WriterField("satellite", (), object, "%5s", 5, "SAT", " "),
        WriterField("used_iode", (), float, "%6d", 6, "IODE", " "),
        WriterField("trans_time_gpsweek", (), object, "%15s", 15, "TRANS_TIME", "wwwwd:ssssss"),
        WriterField("toe_gpsweek", (), object, "%15s", 15, "TOE", "wwwwd:ssssss"),
        WriterField("diff_trans_toe", (), float, "%8d", 8, "TM-TOE"),
        WriterField("age_of_ephemeris", (), float, "%8d", 8, "T-TOE"),
        # WriterField("diff_time_trans",      (),             float,  "%8d",     8, "T-TM",  ),
        WriterField("clk_diff", (), float, "%16.4f", 16, "ΔCLOCK"),
        WriterField("clk_diff_with_dt_mean", (), float, "%16.4f", 16, "ΔCLOCK_MEAN"),
        WriterField("dalong_track", (), float, "%16.4f", 16, "ΔALONG_TRACK"),
        WriterField("dcross_track", (), float, "%16.4f", 16, "ΔCROSS_TRACK"),
        WriterField("dradial", (), float, "%16.4f", 16, "ΔRADIAL"),
        WriterField("orb_diff_3d", (), float, "%16.4f", 16, "ORB_DIFF_3D"),
        WriterField("sisre_orb", (), float, "%16.4f", 16, "SISRE_ORB"),
        WriterField("sisre", (), float, "%16.4f", 16, "SISRE"),
    )

    # Add additional fields used by the writer
    dset.add_text("time_date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime])
    dset.add_text(
        "time_gpsweek", val=[f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.time]
    )
    dset.add_text(
        "trans_time_gpsweek",
        val=[
            f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.used_transmission_time
        ],
    )
    dset.add_text(
        "toe_gpsweek",
        val=[f"{t.gps_ws.week:04.0f}{t.gps_ws.day:1.0f}:{t.gps_ws.seconds:06.0f}" for t in dset.used_toe],
    )
    # dset.add_float("diff_time_trans", val=(dset.time.mjd - dset.used_transmission_time.mjd) * Unit.day2second, Unit="second")
    dset.add_float("dalong_track", val=dset.orb_diff.acr.along, unit=dset.unit("orb_diff.acr.along"))
    dset.add_float("dcross_track", val=dset.orb_diff.acr.cross, unit=dset.unit("orb_diff.acr.cross"))
    dset.add_float("dradial", val=dset.orb_diff.acr.radial, unit=dset.unit("orb_diff.acr.radial"))

    ## Add 'detail' fields used by the writer
    # if write_level <= enums.get_value("write_level", "detail"):
    #    fields += (
    #        WriterField("clk_brdc_com", (), float, "%16.4f", 16, "CLK_BRDC"),
    #        WriterField("clk_precise_com", (), float, "%16.4f", 16, "CLK_PRECISE"),
    #        WriterField("bias_brdc", (), float, "%10.4f", 10, "B_BRDC"),
    #        WriterField("bias_precise", (), float, "%10.4f", 10, "B_PREC"),
    #        WriterField("dt_mean", (), float, "%10.4f", 10, "dt_MEAN"),
    #    )
    #
    #    dset.add_float("dt_mean", val=dset.clk_diff - dset.clk_diff_with_dt_mean, unit="meter")

    # List epochs ordered by satellites
    idx = np.concatenate([np.where(dset.filter(satellite=s))[0] for s in dset.unique("satellite")])

    # Put together fields in an array as specified by the fields-tuple
    output_list = list(zip(*(_get_field(dset, f.field, f.attrs) for f in fields)))
    output_array = np.array(output_list, dtype=[(f.field, f.dtype) for f in fields])[idx]

    # Write to disk
    # NOTE: np.savetxt is used instead of having a loop over all observation epochs, because the performance is better.
    file_path = config.files.path(f"output_sisre_{dset.vars['label']}", file_vars=dset.vars)
    header = [
        _get_header(dset),
        "".join(f"{f.header:>{f.width}s}" for f in fields),
        "".join(f"{f.unit if f.unit else dset.unit(f.field)[0]:>{f.width}s}" for f in fields),
        "_" * sum([f.width for f in fields]),
    ]
    np.savetxt(
        file_path,
        output_array,
        fmt=tuple(f.format for f in fields),
        header="\n".join(header),
        delimiter="",
        encoding="utf8",
    )

    # Append SISRE output path to SISRE output buffer file
    if config.tech.sisre_writer.write_buffer_file.bool:
        sisre_output_buffer.sisre_output_buffer(dset)


def _get_field(dset: "Dataset", field: "str", attrs: Tuple[str]) -> np.ndarray:
    """Get field values of a Dataset specified by the field attributes

    Args:
        dset:     Dataset, a dataset containing the data.
        field:    Field name.
        attrs:    Field attributes (e.g. for Time object: (<scale>, <time format>)).

    Returns:
        Array with Dataset field values
    """
    f = dset[field]
    for attr in attrs:
        f = getattr(f, attr)
    return f


def _get_header(dset: "Dataset") -> str:
    """Get header

    Args:
        dset:   A dataset containing the data.

    Returns:
        Header lines
    """

    pgm = "where " + where.__version__ + "/midgard " + midgard.__version__
    run_by = util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else ""
    file_created = datetime.utcnow().strftime("%Y%m%d %H%M%S") + " UTC"
    header = "PGM: {:s}  RUN_BY: {:s}  DATE: {:s}\n\n".format(pgm, run_by, file_created)
    header += "SISRE ANALYSIS CONFIGURATION\n\n"

    # SISRE configuration
    header += str(config.tech.as_str(key_width=25, width=70, only_used=True)) + "\n\n\n"
    header += _get_paths()

    # Information about used biases and phase center offsets (PCOs)
    header += "{:^4s}{:>10s}{:>10s}{:^26s}{:^26s}\n".format("SAT", "BIAS_BRDC", "BIAS_PREC", "PCO_BRDC", "PCO_PREC")
    header += "{:^4s}{:^10s}{:^10s}{:^26s}{:^26s}\n".format("", "[m]", "[m]", "[m]", "[m]")
    for sat in dset.unique("satellite"):
        header += "{:>3s}{:>10.4f}{:>10.4f}{:>10.4f}{:>8.4f}{:>8.4f}{:>10.4f}{:>8.4f}{:>8.4f}\n".format(
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

    return header + "\n\n"


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
