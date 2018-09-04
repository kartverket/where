"""Write SISRE analysis results

Description:
------------
TODO: Has IGS monitoring group defined a SISRE format?


"""
# Standard library imports
from collections import namedtuple
from datetime import datetime

# External library imports
import midgard
import numpy as np

# Where imports
import where
from where.lib import config
from where.lib import enums
from where.lib import files
from where.lib import plugins
from where.lib.unit import unit
from where.lib import util
from where.writers import sisre_output_buffer

WriterField = namedtuple("WriterField", ["field", "attrs", "dtype", "format", "width", "header", "unit"])
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        field (str):             Dataset field name
        dtype (Numpy dtype):     Type of field
        format (str):            Format string
    """


@plugins.register
def sisre_writer(dset):
    """Write SISRE analysis results

    Args:
        dset:       Dataset, a dataset containing the data.
    """
    write_level = config.tech.get("write_level", default="operational").as_enum("write_level")

    fields = (
        WriterField("time_date", (), object, "%21s", 19, "EPOCH"),
        WriterField("time", ("gps", "mjd"), float, "%14.6f", 14, "", "mjd"),
        WriterField("time_gpsweek", (), object, "%15s", 15, ""),
        WriterField("satellite", (), object, "%5s", 5, "SAT"),
        WriterField("used_iode", (), float, "%6d", 6, "IODE"),
        WriterField("trans_time_gpsweek", (), object, "%15s", 15, "TRANS_TIME"),
        WriterField("toe_gpsweek", (), object, "%15s", 15, "TOE"),
        WriterField("diff_trans_toe", (), float, "%8d", 8, "TM-TOE"),
        WriterField("diff_time_toe", (), float, "%8d", 8, "T-TOE"),
        WriterField("clk_diff_no_mean", (), float, "%16.4f", 16, "ΔCLOCK_NO_MEAN"),
        WriterField("clk_diff", (), float, "%16.4f", 16, "ΔCLOCK"),
        WriterField("dalong_track", (), float, "%16.4f", 16, "ΔALONG_TRACK"),
        WriterField("dcross_track", (), float, "%16.4f", 16, "ΔCROSS_TRACK"),
        WriterField("dradial", (), float, "%16.4f", 16, "ΔRADIAL"),
        WriterField("orb_diff_3d", (), float, "%16.4f", 16, "ORB_DIFF_3D"),
        WriterField("sisre_orb_no_mean", (), float, "%16.4f", 16, "SISRE_ORB_NOMEAN"),
        WriterField("sisre_orb", (), float, "%16.4f", 16, "SISRE_ORB"),
        WriterField("sisre", (), float, "%16.4f", 16, "SISRE"),
    )

    # Add additional fields used by the writer
    dset.add_text(
        "time_date", val=[d.strftime("%Y/%m/%d %H:%M:%S") for d in dset.time.datetime], unit="YYYY/MM/DD hh:mm:ss"
    )
    dset.add_text(
        "time_gpsweek",
        val=[f"{t.gpsweek:04.0f}{t.gpsday:1.0f}:{t.gpssec:06.0f}" for t in dset.time],
        unit="wwwwd:ssssss",
    )
    dset.add_text(
        "trans_time_gpsweek",
        val=[f"{t.gpsweek:04.0f}{t.gpsday:1.0f}:{t.gpssec:06.0f}" for t in dset.used_transmission_time],
        unit="wwwwd:ssssss",
    )
    dset.add_text(
        "toe_gpsweek",
        val=[f"{t.gpsweek:04.0f}{t.gpsday:1.0f}:{t.gpssec:06.0f}" for t in dset.used_toe],
        unit="wwwwd:ssssss",
    )
    dset.add_float(
        "diff_trans_toe", val=(dset.used_transmission_time.mjd - dset.used_toe.mjd) * unit.day2second, unit="second"
    )
    dset.add_float("diff_time_toe", val=(dset.time.mjd - dset.used_toe.mjd) * unit.day2second, unit="second")
    dset.add_float("dalong_track", val=dset.orb_diff_acr.itrs[:, 0], unit=dset.unit("orb_diff_acr.itrs"))
    dset.add_float("dcross_track", val=dset.orb_diff_acr.itrs[:, 1], unit=dset.unit("orb_diff_acr.itrs"))
    dset.add_float("dradial", val=dset.orb_diff_acr.itrs[:, 2], unit=dset.unit("orb_diff_acr.itrs"))

    # Add 'detail' fields used by the writer
    if write_level <= enums.get_value("write_level", "detail"):
        fields += (
            WriterField("clk_brdc_com", (), float, "%16.4f", 16, "CLK_BRDC"),
            WriterField("clk_precise_com", (), float, "%16.4f", 16, "CLK_PRECISE"),
            WriterField("bias_brdc", (), float, "%10.4f", 10, "B_BRDC"),
            WriterField("bias_precise", (), float, "%10.4f", 10, "B_PREC"),
            WriterField("clk_sys", (), float, "%10.4f", 10, "CLK_SYS"),
        )

        dset.add_float("clk_sys", val=dset.clk_diff_no_mean - dset.clk_diff, unit="meter")

    # List epochs ordered by satellites
    idx = np.concatenate([np.where(dset.filter(satellite=s))[0] for s in dset.unique("satellite")])

    # Put together fields in an array as specified by the fields-tuple
    output_list = list(zip(*(_get_field(dset, f.field, f.attrs) for f in fields)))
    output_array = np.array(output_list, dtype=[(f.field, f.dtype) for f in fields])[idx]

    # Write to disk
    # NOTE: np.savetxt is used instead of having a loop over all observation epochs, because the performance is better.
    file_path = files.path(f"output_sisre_{dset.dataset_id}", file_vars=dset.vars)
    header = [
        _get_header(dset),
        "".join(f"{f.header:>{f.width}s}" for f in fields),
        "".join(f"{f.unit if f.unit else dset.unit(f.field):>{f.width}s}" for f in fields),
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


def _get_field(dset, field, attrs):
    """Get field values of a Dataset specified by the field attributes

    Args:
        dset (Dataset):     Dataset, a dataset containing the data.
        field (str):        Field name.
        attrs (tuple):      Field attributes (e.g. for Time object: (<scale>, <time format>)).

    Returns:
        numpy.ndarray:      Array with Dataset field values
    """
    f = dset[field]
    for attr in attrs:
        f = getattr(f, attr)
    return f


def _get_header(dset):
    """Get header

    Args:
        dset (Dataset):       Dataset, a dataset containing the data.

    Returns:
        str:                  Header lines
    """

    pgm = "where " + where.__version__ + "/midgard " + midgard.__version__
    run_by = util.get_user_info()["inst_abbreviation"] if "inst_abbreviation" in util.get_user_info() else ""
    file_created = datetime.utcnow().strftime("%Y%m%d %H%M%S") + " UTC"
    header = ("PGM: {:s}  RUN_BY: {:s}  DATE: {:s}\n\n".format(pgm, run_by, file_created))
    header += "SISRE ANALYSIS RESULTS\n\n"

    # Used input files
    header += str(config.tech.as_str(key_width=25, width=70, only_used=True)) + "\n\n\n"
    header += "{:24s} {}\n".format(
        "Broadcast orbit:", str(files.path("gnss_rinex_nav_M", file_vars=dset.vars))
    )  # TODO: Not correct if station specific file is used!!!
    header += "{:24s} {}\n".format("Precise orbit:", str(files.path("gnss_orbit_sp3", file_vars=dset.vars)))
    if config.tech.clock_product.str == "clk":
        header += "{:24s} {}\n".format("Precise clock:", str(files.path("gnss_rinex_clk", file_vars=dset.vars)))
    header += "{:24s} {}\n".format("Bias:", str(files.path("gnss_sinex_bias", file_vars=dset.vars)))
    header += "{:24s} {}\n".format("ANTEX (broadcast):", str(files.path("sisre_antex_brdc", file_vars=dset.vars)))
    header += "{:24s} {}\n\n\n".format("ANTEX (precise):", str(files.path("gnss_antex", file_vars=dset.vars)))

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
