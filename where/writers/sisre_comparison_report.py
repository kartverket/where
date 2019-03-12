"""Compare different SISRE Where datasets

Description:
------------

A dictionary with datasets is used as input for this writer. The keys of the dictionary have to include the signal type
used for generation of the SISRE dataset. -> TODO: signal type of SISRE dataset has to be defined via meta variable!!!!


Example:
--------

    from where import data
    from where import writers

    # Read a dataset
    dset = data.Dataset(rundate=rundate, tech=tech, stage=stage, dataset_name=name, dataset_id=dataset_id)

    # Write dataset
    writers.write_one('sisre_comparison_report', dset=dsets, do_report=False)

"""
# Standard library imports
from collections import namedtuple
from datetime import datetime
import os
import re
import textwrap

# External library imports
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Midgard imports
from midgard.dev import console
from midgard.dev import plugins

# WHERE imports
import where
from where import apriori
from where import cleaners
from where.lib import config
from where.lib import enums
from where.lib import files
from where.lib import log
from where.lib import gnss
from where.lib import util

FIGURE_FORMAT = "png"
FILE_NAME = __name__.split(".")[-1]


@plugins.register
def sisre_comparison_report(dset):
    """Compare SISRE datasets

    Args:
        dset (list):       List with different SISRE datasets. The datasets contain the data.
    """
    dsets = dset
    df_merged = pd.DataFrame()

    for name, dset in dsets.items():

        if dset.num_obs == 0:
            log.warn(f"Dataset '{name}' is empty.")
            continue

        signal_type = _get_signal_type(dset.meta)
        df = dset.as_dataframe(fields=["satellite", "system", "sisre", "time.gps"])  # , index="time.gps")
        df = df.rename(columns={"sisre": signal_type})

        if df_merged.empty:
            df_merged = df
            continue
        df_merged = df_merged.merge(df, on=["satellite", "system", "time.gps"], how="outer")

    if df_merged.empty:
        log.fatal(f"All given datasets are empty [{', '.join(dsets.keys())}].")

    with files.open(
        file_key="output_sisre_comparison_report", file_vars=dsets[next(iter(dsets))].vars, create_dirs=True, mode="wt"
    ) as fid:
        _header(fid)
        fid.write("#Comparison of SISE analyses\n")
        fid.write("In the following SISE analyses results are compared for:\n\n")
        fid.write("* Monthly 95th percentile SISE for satellites\n")
        fid.write("* Monthly 95th percentile and RMS SISE for signal combinations (users)\n")
        fid.write("\\newpage\n")

        # Generate figure directory to save figures generated for SISRE report
        figure_dir = files.path("output_sisre_comparison_report_figure", file_vars=dset.vars)
        figure_dir.mkdir(parents=True, exist_ok=True)

        fid.write(f"\n\n##Monthly 95th percentile SISE for satellites\n")
        # Produce plot with same yrange than for _plot_bar_sisre_signal_combination_percentile threshold plot
        _plot_bar_sisre_satellite_percentile(
            df_merged, fid, figure_dir, threshold=False, write_table=True, yrange=[0, 2]
        )
        _plot_bar_sisre_satellite_percentile(df_merged, fid, figure_dir, threshold=True, write_table=False)

        fid.write(f"\n\n##Monthly 95th percentile and RMS SISE for signal combinations (users)\n")
        _plot_bar_sisre_signal_combination_percentile(df_merged, fid, figure_dir, threshold=False, write_table=True)
        _plot_bar_sisre_signal_combination_percentile(df_merged, fid, figure_dir, threshold=True, write_table=False)
        _plot_bar_sisre_signal_combination_rms(df_merged, fid, figure_dir, write_table=True)

    # Generate PDF from Markdown file
    _markdown_to_pdf(dset)


def _get_signal_type(meta):
    """Get signal type used for SISRE dataset

    Args:
        meta (dict): Dataset meta dictionary

    Returns:
        str:    Signal type (e.g. E1, E1/E5a, ...)
    """
    if len(meta["systems"]) > 1:
        log.fatal(
            f"The writer '{FILE_NAME}' can only be used, if the dataset is based on one GNSS "
            f"(not '{', '.join(meta['systems'])})'."
        )

    system = meta["systems"][0]
    try:
        signal_type = meta["frequencies"][system].replace("_", "/")

    except KeyError:
        log.fatal(f"No frequencies are defined for GNSS {system!r} for option 'frequencies'.")

    return signal_type


def _header(fid):
    title = "Where SISE comparison"
    user_info = util.get_user_info()
    user = user_info.get("name", user_info["user"])

    fid.write(
        console.dedent(
            f"""---
                title: {title}
                author: Where v{where.__version__} [{user}]
                date: {datetime.now():%Y-%m-%d}
                ---
            """
        )
    )
    fid.write("\\newpage\n\n")


def _markdown_to_pdf(dset):
    """Convert markdown SISRE report file to pdf format

    Args:
       dset (Dataset):           A dataset containing the data.
    """

    if config.where.sisre_report.get("markdown_to_pdf", default=False).bool:
        md_path = str(files.path("output_sisre_comparison_report", file_vars=dset.vars))
        pdf_path = md_path.replace(".md", ".pdf")
        program = "pandoc"

        # Convert markdown to pdf with pandoc
        pandoc_args = ["-f markdown", "-V classoption:twoside", "-N", "-o " + pdf_path, md_path]

        log.info(f"Start: {program} {' '.join(pandoc_args)}")
        status = os.system(f"{program} {' '.join(pandoc_args)}")
        if status != 0:
            log.error(f"{program} failed with error code {status} ({' '.join([program] + pandoc_args)})")

        # TODO: pandoc subprocess call does not work. Why?
        # process = subprocess.run([program] + pandoc_args)
        # if process.returncode:
        #    log.error(f"{program} failed with error code {process.returncode} ({' '.join(process.args)})")


def _plot_bar_sisre_signal_combination_percentile(df, fid, figure_dir, threshold=False, write_table=False):
    """Generate bar plot with monthly SISRE 95% percentile for each GNSS signal combination

    Args:
       df (Dataframe):          Dataframe with time, satellite, system and GNSS signal combinations as columns
       fid (_io.TextIOWrapper): File object.
       figure_dir (PosixPath):  Figure directory.
       threshold (bool):        Plot threshold.
       write_table (bool):      Write table.
    """
    # df_monthly_percentile = df.set_index('time.gps').resample('D', how=lambda x: np.nanpercentile(x, q=95))
    # df_monthly_percentile = df.dropna().set_index("time.gps").resample("M", how=lambda x: np.percentile(x, q=95))
    df_monthly_percentile = df.set_index("time.gps").resample("M", how=lambda x: np.nanpercentile(x, q=95))
    df_monthly_percentile.index = df_monthly_percentile.index.strftime("%b-%Y")
    df_monthly_percentile.transpose().plot(kind="bar")

    if write_table:
        _write_dataframe_to_markdown(
            fid,
            df_monthly_percentile.transpose(),
            float_format="6.3f",
            caption=f"95th percentile SISE results for signal combinations in meter: ",
        )

    if threshold:
        plt.axhline(2, color="r")
        plot_filename = f"plot_bar_sisre_signal_combination_percentile_threshold.{FIGURE_FORMAT}"
    else:
        plot_filename = f"plot_bar_sisre_signal_combination_percentile.{FIGURE_FORMAT}"
    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISE [m] (95th percentile)")
    # plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1) #legend on right side
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3) #legend on top
    plt.legend(bbox_to_anchor=(0.8, -0.15), loc=1, borderaxespad=0.0, ncol=3)  # legend below
    # plt.legend(bbox_to_anchor=(0.8, -0.1), loc=1, borderaxespad=0., ncol=3) #legend below without xlabel
    plt.tight_layout()
    plt.savefig(figure_dir / plot_filename)
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly 95th percentile of global average SISE for single- and dual-frequency users]({figure_dir}/{plot_filename})\n"
    )
    fid.write("\\newpage\n")


def _plot_bar_sisre_signal_combination_rms(df, fid, figure_dir, write_table=False):
    """Generate bar plot with monthly SISRE RMS for each GNSS signal combination

    Args:
       df (Dataframe):          Dataframe with time, satellite, system and GNSS signal combinations as columns
       fid (_io.TextIOWrapper): File object.
       figure_dir (PosixPath):  Figure directory.
       write_table (bool):      Write table.
    """
    df_monthly_rms = df.set_index("time.gps").resample("M", how=lambda x: np.sqrt(np.nanmean(np.square(x))))
    df_monthly_rms.index = df_monthly_rms.index.strftime("%b-%Y")
    df_monthly_rms.transpose().plot(kind="bar")

    if write_table:
        _write_dataframe_to_markdown(
            fid,
            df_monthly_rms.transpose(),
            float_format="6.3f",
            caption=f"RMS SISE results for signal combinations in meter: ",
        )

    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISE RMS [m]")
    # plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1) #legend on right side
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3) #legend on top
    plt.legend(bbox_to_anchor=(0.8, -0.15), loc=1, borderaxespad=0.0, ncol=3)  # legend below
    # plt.legend(bbox_to_anchor=(0.8, -0.1), loc=1, borderaxespad=0., ncol=3) #legend below without xlabel
    plt.tight_layout()
    plt.savefig(figure_dir / f"plot_bar_sisre_signal_combination_rms.{FIGURE_FORMAT}")
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly RMS of global average SISE for single- and dual-frequency users]({figure_dir}/plot_bar_sisre_signal_combination_rms.{FIGURE_FORMAT})\n"
    )
    fid.write("\\newpage\n")


def _plot_bar_sisre_satellite_percentile(df, fid, figure_dir, threshold=False, write_table=False, yrange=None):
    """Generate bar plot with monthly SISRE 95% percentile for each satellite

    Args:
       df (Dataframe):          Dataframe with time, satellite, system and GNSS signal combinations as columns
       fid (_io.TextIOWrapper): File object.
       figure_dir (PosixPath):  Figure directory.
       write_table (bool):      Write table.
    """
    # Get user types by keeping order
    user_types = list(df.columns)
    for value in ["satellite", "system", "time.gps"]:
        user_types.remove(value)

    fig, axes = plt.subplots(len(user_types), 1, sharex=False, sharey=True)  # figsize=(6, 6));

    # Guarantee that 'axes' is iterable, which is needed by using 'zip' command.
    # Note: 'plt.subplots' does not return an iterable numpy array, if row and column number is one. This is the case,
    #       if 'user_types' has only one element.
    if not isinstance(axes, (np.ndarray)):
        axes = np.array([axes])

    if threshold:
        plot_filename = f"plot_bar_sisre_satellite_percentile_threshold.{FIGURE_FORMAT}"
    else:
        plot_filename = f"plot_bar_sisre_satellite_percentile.{FIGURE_FORMAT}"

    # plt.subplots_adjust(wspace=0.5, hspace=0.5);
    fig.set_figheight(10)  # inches

    # Generate subplots
    for ax, user in zip(axes, user_types):
        df_user = df.pivot(index="time.gps", columns="satellite", values=user)
        df_user_monthly_percentile = df_user.resample("M", how=lambda x: np.nanpercentile(x, q=95))
        df_user_monthly_percentile.index = df_user_monthly_percentile.index.strftime("%b-%Y")
        df_user_monthly_percentile.transpose().plot(kind="bar", ax=ax, legend=False, title=user)
        ax.set_ylabel("SISE [m] (95th percentile)")
        if threshold:
            ax.axhline(7, color="r")
        if yrange is not None:
            ax.set_ylim(yrange)

        if write_table:
            _write_dataframe_to_markdown(
                fid,
                df_user_monthly_percentile.transpose(),
                float_format="6.3f",
                caption=f"95th percentile SISE results for signal combination **{user}** in meter: ",
            )
            fid.write("\\newpage\n")

    plt.legend(bbox_to_anchor=(0.8, -0.3), loc=1, borderaxespad=0.0, ncol=3)
    plt.tight_layout()
    plt.savefig(figure_dir / plot_filename)
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly 95th percentile of global average SISE of each satellite for single- and dual-frequency users]({figure_dir}/{plot_filename})\n"
    )
    fid.write("\\newpage\n")


def _write_dataframe_to_markdown(fid, df, float_format="", caption=""):
    """Write Pandas DataFrame to Markdown

    Args:
        fid (_io.TextIOWrapper):    File object
        df (DataFrame):             Pandas DataFrame
        float_format (str):         Define formatters for float columns
    """
    column_length = [len(c) for c in df.columns]

    # Write caption
    fid.write(caption + "\n")

    # Write header
    if list(df.index):  # Add DataFrame index to header
        num_space = len(max(df.index))
        head_line_1 = "\n| {} ".format(" " * num_space)
        head_line_2 = "|-{}-".format("-" * num_space)
    else:
        header_1 = ""

    fid.write(head_line_1 + "| {} |\n".format(" | ".join(list(df.columns))))
    fid.write(head_line_2 + "|-{}:|\n".format("-|-".join([n * "-" for n in column_length])))

    # Write data
    for index, row in df.iterrows():

        line = "| {idx:s} |".format(idx=index) if index else ""  # Add DataFrame index column

        for _, v in row.items():
            if isinstance(v, float):
                line = line + " {:{fmt}} |".format(v, fmt=float_format)
            else:
                line = line + " {} |".format(v)
        fid.write(line + "\n")
    fid.write("\n")
