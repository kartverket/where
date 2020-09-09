"""Compare two different GNSS Where datasets

Description:
------------

A dictionary with datasets is used as input for this writer. A differenced dataset 'dset_diff' is generated, which
includes difference between 1st and 2nd dataset for common fields (except Time object fields, which are excluded). In
addition plots and a report is generated.

Example:
--------
where_tools compare 2019 7 1 2019 7 1 --gnss --writers=gnss_compare_datasets --items='estimate, springarp' --specifier='stage' --dset_name=krss -T


"""
# Standard library imports
from collections import namedtuple
from typing import Any, Dict, List, Set

# External library imports
import matplotlib.pyplot as plt
import numpy as np

# Midgard imports
from midgard.data._time import GpsTime, UtcTime
from midgard.dev import plugins
from midgard.plot.matplotlib_extension import plot_scatter_subplots
from midgard.math.unit import Unit

# Where imports
import where
from where.lib import config
from where.lib import log
from where.writers._report import Report

FIGURE_FORMAT = "png"
FILE_NAME = __name__.split(".")[-1]


WriterField = namedtuple("WriterField", ["label", "title", "unit"])
WriterField.__new__.__defaults__ = (None,) * len(WriterField._fields)
WriterField.__doc__ = """A convenience class for defining a output field for the writer

    Args:
        label (str):     Y-axis label name of plot
        title (str):     Plot title
        unit  (str):     Unit of field           
    """

FIELDS = {
    "gnss_ionosphere": WriterField("Delay", "Ionosphere delay", "m"),
    "gnss_range": WriterField("Range", "Range between satellite and receiver", "m"),
    "gnss_satellite_clock": WriterField("Clock", "Satellite clock correction", "m"),
    "gnss_total_group_delay": WriterField("TGD", "Total group delay", "m"),
    "gdop": WriterField("GDOP", "Geometric dilution of precision", "m"),
    "hdop": WriterField("HDOP", "Horizontal dilution of precision", "m"),
    "num_satellite_available": WriterField("#satellites", "Number of available satellites", ""),
    "num_satellite_used": WriterField("#satellites", "Number of used satellites", ""),
    "sat_pos_x": WriterField("X", "Satellite position", "m"),
    "sat_pos_y": WriterField("Y", "Satellite position", "m"),
    "sat_pos_z": WriterField("Z", "Satellite position", "m"),
    "residual": WriterField("Residual", "Residual", "m/s"),
    "site_pos_vs_ref_east": WriterField("East", "Site position vs. reference position", "m"),
    "site_pos_vs_ref_north": WriterField("North", "Site position vs. reference position", "m"),
    "site_pos_vs_ref_up": WriterField("Up", "Site position vs. reference position", "m"),
    "site_vel_3d": WriterField("3D", "3D site velocity", "m/s"),
    "site_vel_h": WriterField("HV", "Horizontal site velocity", "m/s"),
    "site_vel_v": WriterField("VV", "Vertical site velocity", "m/s"),
    "site_vel_x": WriterField("X", "Site velocity", "m/s"),
    "site_vel_y": WriterField("Y", "Site velocity", "m/s"),
    "site_vel_z": WriterField("Z", "Site velocity", "m/s"),
    "site_vel_east": WriterField("East", "Site velocity", "m/s"),
    "site_vel_north": WriterField("North", "Site velocity", "m/s"),
    "site_vel_up": WriterField("Up", "Site velocity", "m/s"),
    "tdop": WriterField("TDOP", "Time dilution of precision", "m"),
    "troposphere_dT": WriterField("Delay", "Tropospheric delay", "m"),
    "vdop": WriterField("VDOP", "Vertical dilution of precision", "m"),
}


@plugins.register
def gnss_compare_datasets(dset: Dict[str, "Dataset"]) -> None:
    """Compare two different GNSS Where datasets

    Args:
        dset:  Dictionary with station name as keys and the belonging Dataset as value
    """

    dset1 = dset[list(dset.keys())[0]]
    dset2 = dset[list(dset.keys())[1]]

    # Decimate datasets
    difference_by = _get_difference_by(dset1.fields, dset2.fields)
    _decimate_datasets(dset1, dset2, difference_by)
    # +MURKS
    vars_ = dset2.vars.copy()
    vars_["stage"] = vars_["stage"] + "x"
    dset2_path = config.files.path("dataset", file_vars=vars_)
    if dset2_path.exists():
        dset2_path.unlink()
    dset2.write_as(stage=vars_["stage"])
    # -MURKS
    if dset1.num_obs == 0 or dset2.num_obs == 0:
        log.fatal(
            f"Nothing to compare. Number of observations are zero at least for one dataset "
            f"(dset1: {dset1.num_obs}, dset2: {dset2.num_obs})."
        )

    # Get common fields in both Datasets
    common_fields = _get_common_fields(dset1, dset2)

    # Generate difference of datasets
    ddiff = _difference_datasets(dset1, dset2, difference_by)

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path("output_gnss_vel_report_figure", file_vars=dset1.vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _plot(dset1, dset2, ddiff, common_fields, figure_dir)

    # Generate report
    path = config.files.path("output_gnss_vel_report", file_vars=dset1.vars)
    with config.files.open_path(path, create_dirs=True, mode="wt") as fid:
        rpt = Report(fid, rundate=dset1.vars["rundate"], path=path, description="Comparison of two GNSS datasets")
        rpt.title_page()
        _add_to_report(rpt, common_fields, figure_dir)
        rpt.markdown_to_pdf()


def _add_to_report(rpt: "Report", common_fields: Set[str], figure_dir: "pathlib.PosixPath") -> None:
    """Add figures to report

    Args:
        rpt:            Report object.
        common_fields:  Set with fields common in both datasets
        figure_dir:     Figure directory.
    """
    dset1_name = config.where.gnss_compare_datasets.get("dset1_name", default="dset1").str
    dset2_name = config.where.gnss_compare_datasets.get("dset2_name", default="dset2").str

    for field in common_fields:
        rpt.add_figure(
            f"{figure_dir}/plot_{field}.{FIGURE_FORMAT}",
            caption=f"Difference between {dset1_name} and {dset2_name} dataset field **{field}**.",
            clearpage=True,
        )


def _plot(
    dset1: "Dataset", dset2: "Dataset", ddiff: "Dataset", common_fields: Set[str], figure_dir: "pathlib.PosixPath"
) -> None:
    """Generate plots

    Args:
        dset1:          First dataset containing the data.
        dset2:          Second dataset containing the data.
        ddiff:          Dataset containing the differences for each field between Dataset 'dset1' and Dataset 'dset2'
        common_fields:  Set with fields common in both datasets
        figure_dir:     Figure directory.
    """
    dset1_name = config.where.gnss_compare_datasets.get("dset1_name", default="dset1").str + ":"
    dset2_name = config.where.gnss_compare_datasets.get("dset2_name", default="dset2").str + ":"

    for field in common_fields:
        ylabel = FIELDS[field].label if field in FIELDS.keys() else field.lower()
        title = FIELDS[field].title if field in FIELDS.keys() else ""
        unit = Unit(FIELDS[field].unit).units if field in FIELDS.keys() else Unit(dset1.unit(field)).units
        options = _set_plot_config(title=title)

        plot_scatter_subplots(
            x_array=dset1.time.gps.datetime,
            y_arrays=[dset1[field], dset2[field], ddiff[field]],
            xlabel="Time [GPS]",
            ylabels=[f"{dset1_name} {ylabel}", f"{dset2_name} {ylabel}", f"Difference: {ylabel}"],
            colors=["steelblue", "darkorange", "limegreen"],
            y_units=[f"{unit:~P}", f"{unit:~P}", f"{unit:~P}"],
            figure_path=figure_dir / f"plot_{field}.{FIGURE_FORMAT}",
            opt_args=options,
        )


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
    import IPython; IPython.embed()
    for values in dset.values(*fields):
        line = ""
        for value in values:
            if isinstance(value, (GpsTime, UtcTime)):  # TODO: Check if more time datatypes should be defined.
                line = f"{line}{value.gps.isot}"
            else:
                line = f"{line}{value}"

        concatenated_fields.append(line.strip())

    return concatenated_fields


def _decimate_datasets(dset1: "Dataset", dset2: "Dataset", decimate_by: List[str]) -> None:
    """Decimate given datasets, that they only incluce corresponding observations

    Args:
        dset1:        First dataset containing the data.
        dset2:        Second dataset containing the data.
        decimate_by:  List with field name to use for decimating given datasets
    """
    keep_idx1 = np.zeros(dset1.num_obs, dtype=bool)
    keep_idx2 = np.zeros(dset2.num_obs, dtype=bool)

    #+TODO
    # Change time scale from GPS to UTC
    dset1_utc = dset1.time.utc
    del dset1.time
    dset1.add_time("time", val= dset1_utc, scale="utc", fmt="datetime")

    dset2_utc = dset2.time.utc
    del dset2.time
    dset2.add_time("time", val= dset2_utc, scale="utc", fmt="datetime")
    #-TODO

    # Get common dataset data   
    dset1_index_data = [dset1[n.strip()] for n in decimate_by]
    dset2_index_data = [dset2[n.strip()] for n in decimate_by]
    A = np.rec.fromarrays(dset1_index_data)
    B = np.rec.fromarrays(dset2_index_data)
    common, dset1_idx, dset2_idx = np.intersect1d(A, B, return_indices=True)
    
    keep_idx1[dset1_idx] = True
    keep_idx2[dset2_idx] = True
    
    # Decimate datasets
    dset1.subset(keep_idx1)
    dset2.subset(keep_idx2)


def _difference_datasets(dset1: "Dataset", dset2: "Dataset", difference_by: List[str]) -> "Dataset":
    """Generate difference between given datasets by using defined fields

    Args:
        dset1:          First dataset containing the data.
        dset2:          Second dataset containing the data.
        difference_by:  List with field name to use for differencing given datasets

    Returns:
        ddiff: Dataset containing the differences for each common field between Dataset 'dset1' and Dataset 'dset2'
    """
    ddiff = dset1.difference(dset2, index_by=",".join(difference_by))
    ddiff.vars.update(dset1.vars)
    ddiff.write_as(pipeline=dset1.vars["pipeline"], stage="dset_diff")
    return ddiff


def _get_common_fields(dset1: "Dataset", dset2: "Dataset") -> Set[str]:
    """Get list with common fields in both Datasets and prepare datasets for plotting

    Args:
        dset1:   First dataset containing the data.
        dset2:   Second dataset containing the data.
    
    Returns:
        Set with common fields
    """
    common_fields = set(dset1.fields) & set(dset2.fields)
    common_fields = (
        set(common_fields) & set(config.where.gnss_plot.fields.list)
        if config.where.gnss_plot.fields.list
        else common_fields
    )

    # Remove Time object and text fields
    remove_fields = set()
    for field in common_fields:
        if "." in field:
            remove_fields.add(field)
            continue
        elif isinstance(dset1[field], (GpsTime, UtcTime)):  # TODO: Check if more time datatypes should be defined.
            remove_fields.add(field)
            continue
        elif field in ["satellite", "system"]:
            remove_fields.add(field)
            continue

    # Add satellite position coordinate fields to dataset
    if "sat_pos" in common_fields:
        dset1.add_float("sat_pos_x", val=dset1.sat_pos.itrs[:, 0])
        dset1.add_float("sat_pos_y", val=dset1.sat_pos.itrs[:, 1])
        dset1.add_float("sat_pos_z", val=dset1.sat_pos.itrs[:, 2])
        del dset1.sat_pos

        dset2.add_float("sat_pos_x", val=dset2.sat_pos.itrs[:, 0])
        dset2.add_float("sat_pos_y", val=dset2.sat_pos.itrs[:, 1])
        dset2.add_float("sat_pos_z", val=dset2.sat_pos.itrs[:, 2])
        del dset2.sat_pos

        remove_fields.add("sat_pos")
        common_fields.update(("sat_pos_x", "sat_pos_y", "sat_pos_z"))

    # Remove unnecessary fields
    common_fields -= remove_fields

    if not common_fields:
        log.fatal("Nothing to compare. No common fields in datasets.")

    return common_fields


def _get_difference_by(fields1: List[str], fields2: List[str]) -> List[str]:
    """Get list with common fields used to decimate and difference given Datasets

    Args:
        fields1: Fields of 1st Dataset
        fields2: Fields of 2nd Dataset
    
    Returns:
        List with common fields to decimate and difference given Datasets
    """
    difference_by = []
    common_fields = set(fields1) & set(fields2)

    for field in ["time", "satellite"]:
        if field in common_fields:
            difference_by.append(field)

    return difference_by


def _set_plot_config(title: str) -> Dict[str, Any]:
    """Set matplotlib configuration by reading gnss_plot configuration

    Args:
        title: Title used as default, if not given in configuration file. 

    Returns:
        Matplotlib configuration options
    """
    # Change fontsize of labels
    fontsize = config.where.gnss_plot.get("fontsize", default=10).str
    fontsize = 9 if not fontsize else int(fontsize)
    plt.rcParams["axes.titlesize"] = fontsize
    plt.rcParams["axes.labelsize"] = fontsize
    plt.rcParams["font.size"] = fontsize
    plt.rcParams["figure.titlesize"] = fontsize
    plt.rcParams["legend.fontsize"] = fontsize
    plt.rcParams["xtick.labelsize"] = fontsize
    plt.rcParams["ytick.labelsize"] = fontsize

    # Define additional matplotlib/plot configuration
    options = {
        "alpha": config.where.gnss_plot.get("alpha", default=1).int,
        "dpi": config.where.gnss_plot.get("dpi", default=200).int,
        "color": config.where.gnss_plot.get("color", default="").str,
        "colormap": config.where.gnss_plot.get("colormap", default="tab10").str,
        "figsize": tuple([int(v) for v in config.where.gnss_plot.get("figsize", default="6,4").list]),
        "fontsize": fontsize,
        "legend": config.where.gnss_plot.get("legend", default=True).bool,
        "marker": config.where.gnss_plot.get("marker", default=".").str,
        "markersize": config.where.gnss_plot.get("markersize", default=9).int,
        "plot_to": "file",
        "sharey": False,
        "statistic": ["rms", "mean", "std", "min", "max", "percentile"],
        "subplot": config.where.gnss_plot.get("subplot", default=True).bool,
        "title": config.where.gnss_plot.get("title", default=" ").str
        if config.where.gnss_plot.get("title", default=" ").str
        else title,
    }

    return options
