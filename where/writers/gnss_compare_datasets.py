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
import midgard
from midgard.data._time import GpsTime, UtcTime
from midgard.data import collection
from midgard.dev import plugins
from midgard.plot.matplotext import MatPlotExt
from midgard.math.unit import Unit

# Where imports
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
    "clk_diff_dt_mean": WriterField("avg(clock)", "Epochwise satellite clock average (dB_mean)", "m"),
    "clk_diff": WriterField("Δclock", "Satellite clock difference", "m"),
    "clk_diff_with_dt_mean": WriterField("Δclock (mean)", "Satellite clock difference (dH_mean)", "m"),
    "delay.gnss_ionosphere": WriterField("Delay", "Ionosphere delay", "m"),
    "delay.gnss_range": WriterField("Range", "Range between satellite and receiver", "m"),
    "delay.gnss_satellite_clock": WriterField("Clock", "Satellite clock correction", "m"),
    "delay.gnss_total_group_delay": WriterField("TGD", "Total group delay", "m"),
    "dradial": WriterField("Δradial", "Orbit difference - radial", "m"),
    "gdop": WriterField("GDOP", "Geometric dilution of precision", ""),
    "hdop": WriterField("HDOP", "Horizontal dilution of precision", ""),
    "num_satellite_available": WriterField("#satellites", "Number of available satellites", ""),
    "num_satellite_used": WriterField("#satellites", "Number of used satellites", ""),
    "pdop": WriterField("PDOP", "Position dilution of precision", ""),
    "orb_diff_3d": WriterField("3D orbit error", "3D orbit error", "m"),
    "orb_diff_x": WriterField("X", "Orbit difference - X", "m"),
    "orb_diff_y": WriterField("Y", "Orbit difference - Y", "m"),
    "orb_diff_z": WriterField("Z", "Orbit difference - Z", "m"),
    "residual": WriterField("Residual", "Residual", "m"),
    "sat_pos_x": WriterField("X", "Satellite position", "m"),
    "sat_pos_y": WriterField("Y", "Satellite position", "m"),
    "sat_pos_z": WriterField("Z", "Satellite position", "m"),
    "sat_vel_x": WriterField("VX", "Satellite velocity", "m/s"),
    "sat_vel_y": WriterField("VY", "Satellite velocity", "m/s"),
    "sat_vel_z": WriterField("VZ", "Satellite velocity", "m/s"),
    "sisre": WriterField("SISE", "Signal-in-space ranging error (URE_Av_mean)", "m"),
    "site_pos_vs_ref_east": WriterField("East", "Site position vs. reference position", "m"),
    "site_pos_vs_ref_north": WriterField("North", "Site position vs. reference position", "m"),
    "site_pos_vs_ref_up": WriterField("Up", "Site position vs. reference position", "m"),
    "site_vel_3d": WriterField("3D", "3D site velocity", "m/s"),
    "site_vel_h": WriterField("HV", "Horizontal site velocity", "m/s"),
    "site_vel_v": WriterField("VV", "Vertical site velocity", "m/s"),
    "site_vel_x": WriterField("VX", "Site velocity", "m/s"),
    "site_vel_y": WriterField("VY", "Site velocity", "m/s"),
    "site_vel_z": WriterField("VZ", "Site velocity", "m/s"),
    "site_vel_east": WriterField("East", "Site velocity", "m/s"),
    "site_vel_north": WriterField("North", "Site velocity", "m/s"),
    "site_vel_up": WriterField("Up", "Site velocity", "m/s"),
    "sqrt_a2_c2": WriterField("SQRT(Δa^2 + Δc^2)", "Orbit difference (dAC)", "m"),
    "tdop": WriterField("TDOP", "Time dilution of precision", ""),
    "troposphere_dT": WriterField("Delay", "Tropospheric delay", "m"),
    "used_iode": WriterField("IODE", "Issue of ephemeris data", ""),
    "vdop": WriterField("VDOP", "Vertical dilution of precision", ""),
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
    if dset1.num_obs == 0 or dset2.num_obs == 0:
        log.fatal(
            f"Nothing to compare. Number of observations are zero at least for one dataset "
            f"(dset1: {dset1.num_obs}, dset2: {dset2.num_obs})."
        )
        
    # Get common fields in both Datasets (+ adding of necessary fields)
    common_fields = _get_common_fields(dset1, dset2)

    # Generate difference of datasets
    ddiff = _difference_datasets(dset1, dset2, difference_by)

    # Generate figure directory to save figures generated for GNSS report
    figure_dir = config.files.path(f"output_{dset1.vars['pipeline']}_report_figure", file_vars=dset1.vars)
    figure_dir.mkdir(parents=True, exist_ok=True)

    # Generate plots
    _plot(dset1, dset2, ddiff, common_fields, figure_dir)

    # Generate report
    path = config.files.path(f"output_{dset1.vars['pipeline']}_report", file_vars=dset1.vars)
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

    for field in sorted(common_fields):

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
        if field in FIELDS.keys():
            unit = f"{Unit(FIELDS[field].unit).units:~P}"
        else:
            if dset1.unit(field):
                unit = f"{Unit(dset1.unit(field)[0]).units:~P}"
            else:
                unit = ""
        options = _set_plot_config(title=title)

        plt = MatPlotExt()
        plt.plot_subplots(
            x_array=dset1.time.gps.datetime,
            y_arrays=[dset1[field], dset2[field], ddiff[field]],
            xlabel="Time [GPS]",
            ylabels=[f"{dset1_name} {ylabel}", f"{dset2_name} {ylabel}", f"Difference: {ylabel}"],
            colors=["steelblue", "darkorange", "limegreen"],
            y_units=[unit, unit, unit],
            figure_path=figure_dir / f"plot_{field}.{FIGURE_FORMAT}",
            options=options,
        )


def _decimate_datasets(dset1: "Dataset", dset2: "Dataset", decimate_by: List[str]) -> None:
    """Decimate given datasets, that they only include corresponding observations

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
    # intersect1d does not like to compare unicode strings of different lengths
    # use object as dtype for all fields to avoid the problem
    dtype_dset1 = ",".join("O"*len(dset1_index_data))
    dtype_dset2 = ",".join("O"*len(dset2_index_data))

    A = np.rec.fromarrays(dset1_index_data, dtype=dtype_dset1)
    B = np.rec.fromarrays(dset2_index_data, dtype=dtype_dset2)
        
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
        #if "." in field:
        #    remove_fields.add(field)
        #    continue
        if isinstance(dset1[field], (GpsTime, UtcTime, collection.Collection)):  # TODO: Check if more time datatypes should be defined.
            remove_fields.add(field)
            continue

        # Skip text fields
        if field in dset1._fields.keys():
            if type(dset1._fields[field]) == midgard.data.fieldtypes.text.TextField:
                remove_fields.add(field)
                continue
        
    # Add Position fields to dataset
    if "orb_diff" in common_fields:
        
        dset1.add_float("orb_diff_x", val=dset1.orb_diff.pos.trs.x)
        dset1.add_float("orb_diff_y", val=dset1.orb_diff.pos.trs.y)
        dset1.add_float("orb_diff_z", val=dset1.orb_diff.pos.trs.z)
        del dset1.orb_diff

        dset2.add_float("orb_diff_x", val=dset2.orb_diff.pos.trs.x)
        dset2.add_float("orb_diff_y", val=dset2.orb_diff.pos.trs.y)
        dset2.add_float("orb_diff_z", val=dset2.orb_diff.pos.trs.z)
        del dset2.orb_diff

        remove_fields.add("orb_diff")
        common_fields.update(("orb_diff_x", "orb_diff_y", "orb_diff_z"))

    if "sat_pos" in common_fields:
        dset1.add_float("sat_pos_x", val=dset1.sat_pos.trs.x)
        dset1.add_float("sat_pos_y", val=dset1.sat_pos.trs.y)
        dset1.add_float("sat_pos_z", val=dset1.sat_pos.trs.z)
        del dset1.sat_pos

        dset2.add_float("sat_pos_x", val=dset2.sat_pos.trs.x)
        dset2.add_float("sat_pos_y", val=dset2.sat_pos.trs.y)
        dset2.add_float("sat_pos_z", val=dset2.sat_pos.trs.z)
        del dset2.sat_pos

        remove_fields.add("sat_pos")
        common_fields.update(("sat_pos_x", "sat_pos_y", "sat_pos_z"))

    if "site_pos" in common_fields:
        dset1.add_float("site_pos_x", val=dset1.site_pos.trs.x, unit="meter")
        dset1.add_float("site_pos_y", val=dset1.site_pos.trs.y, unit="meter")
        dset1.add_float("site_pos_z", val=dset1.site_pos.trs.z, unit="meter")
        del dset1.site_pos

        dset2.add_float("site_pos_x", val=dset2.site_pos.trs.x, unit="meter")
        dset2.add_float("site_pos_y", val=dset2.site_pos.trs.y, unit="meter")
        dset2.add_float("site_pos_z", val=dset2.site_pos.trs.z, unit="meter")
        del dset2.site_pos

        remove_fields.add("site_pos")
        common_fields.update(("site_pos_x", "site_pos_y", "site_pos_z"))

    # Add PosVel fields to dataset
    if "gnss_earth_rotation" in common_fields:
        
        dset1.add_float("gnss_earth_rotation_x", val=dset1.gnss_earth_rotation.pos.trs.x)
        dset1.add_float("gnss_earth_rotation_y", val=dset1.gnss_earth_rotation.pos.trs.y)
        dset1.add_float("gnss_earth_rotation_z", val=dset1.gnss_earth_rotation.pos.trs.z)
        dset1.add_float("gnss_earth_rotation_vx", val=dset1.gnss_earth_rotation.vel.trs.x)
        dset1.add_float("gnss_earth_rotation_vy", val=dset1.gnss_earth_rotation.vel.trs.y)
        dset1.add_float("gnss_earth_rotation_vz", val=dset1.gnss_earth_rotation.vel.trs.z)
        del dset1.gnss_earth_rotation

        dset2.add_float("gnss_earth_rotation_x", val=dset2.gnss_earth_rotation.pos.trs.x)
        dset2.add_float("gnss_earth_rotation_y", val=dset2.gnss_earth_rotation.pos.trs.y)
        dset2.add_float("gnss_earth_rotation_z", val=dset2.gnss_earth_rotation.pos.trs.z)
        dset2.add_float("gnss_earth_rotation_vx", val=dset2.gnss_earth_rotation.vel.trs.x)
        dset2.add_float("gnss_earth_rotation_vy", val=dset2.gnss_earth_rotation.vel.trs.y)
        dset2.add_float("gnss_earth_rotation_vz", val=dset2.gnss_earth_rotation.vel.trs.z)
        del dset2.gnss_earth_rotation

        remove_fields.add("gnss_earth_rotation")
        common_fields.update(("gnss_earth_rotation_x", "gnss_earth_rotation_y", "gnss_earth_rotation_z", "gnss_earth_rotation_vx", "gnss_earth_rotation_vy", "gnss_earth_rotation_vz"))

    if "sat_posvel" in common_fields:
        
        dset1.add_float("sat_pos_x", val=dset1.sat_posvel.pos.trs.x)
        dset1.add_float("sat_pos_y", val=dset1.sat_posvel.pos.trs.y)
        dset1.add_float("sat_pos_z", val=dset1.sat_posvel.pos.trs.z)
        dset1.add_float("sat_vel_x", val=dset1.sat_posvel.vel.trs.x)
        dset1.add_float("sat_vel_y", val=dset1.sat_posvel.vel.trs.y)
        dset1.add_float("sat_vel_z", val=dset1.sat_posvel.vel.trs.z)
        del dset1.sat_posvel

        dset2.add_float("sat_pos_x", val=dset2.sat_posvel.pos.trs.x)
        dset2.add_float("sat_pos_y", val=dset2.sat_posvel.pos.trs.y)
        dset2.add_float("sat_pos_z", val=dset2.sat_posvel.pos.trs.z)
        dset2.add_float("sat_vel_x", val=dset2.sat_posvel.vel.trs.x)
        dset2.add_float("sat_vel_y", val=dset2.sat_posvel.vel.trs.y)
        dset2.add_float("sat_vel_z", val=dset2.sat_posvel.vel.trs.z)
        del dset2.sat_posvel

        remove_fields.add("sat_posvel")
        common_fields.update(("sat_pos_x", "sat_pos_y", "sat_pos_z", "sat_vel_x", "sat_vel_y", "sat_vel_z"))
        

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

    for field in ["time", "satellite", "system"]:
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
        "figsize": tuple([int(v) for v in config.where.gnss_plot.get("figsize", default="8,6").list]),
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
