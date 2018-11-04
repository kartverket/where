"""Compare different SISRE Where datasets

Description:
------------

asdf




"""
# Standard library imports
from collections import namedtuple
from datetime import datetime
import getpass
import os
import re
import textwrap

# External library imports
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# WHERE imports
import where
from where import apriori
from where import cleaners
from where.lib import config
from where.lib import enums
from where.lib import files
from where.lib import log
from where.lib import gnss
from where.lib import plugins

FIGURE_FORMAT = "png"


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

        user_type_name = _get_user_type_name(name)
        df = dset.as_dataframe(fields=["satellite", "system", "sisre", "time.gps"])  # , index="time.gps")
        df = df.rename(columns={"sisre": user_type_name})

        if df_merged.empty:
            df_merged = df
            continue
        df_merged = df_merged.merge(df, on=["satellite", "system", "time.gps"], how="outer")

    if df_merged.empty:
        log.fatal(f"All given datasets are empty [{', '.join(dsets.keys())}].")

    with files.open(
        file_key="output_sisre_comparison_report", file_vars=dsets[next(iter(dsets))].vars, mode="wt"
    ) as fid:
        _header(fid)
        fid.write("#Comparison of SISRE analyses\n")

        # Generate figure directory to save figures generated for SISRE report
        figure_dir = files.path("output_sisre_comparison_report_figure", file_vars=dset.vars)
        figure_dir.mkdir(parents=True, exist_ok=True)

        _plot_bar_sisre_satellite_percentile(df_merged, fid, figure_dir, threshold=False)
        _plot_bar_sisre_satellite_percentile(df_merged, fid, figure_dir, threshold=True)
        _plot_bar_sisre_signal_combination_percentile(df_merged, fid, figure_dir, threshold=False)
        _plot_bar_sisre_signal_combination_percentile(df_merged, fid, figure_dir, threshold=True)
        _plot_bar_sisre_signal_combination_rms(df_merged, fid, figure_dir)

    # Generate PDF from Markdown file
    _markdown_to_pdf(dset)


def _get_user_type_name(dset_name):
    # TODO: Pattern search inav_e1 does not work. How to solve?
    user_type_name_def = {"fnav_e1e5a": "E1E5a", "inav_e1_": "E1", "inav_e1e5b": "E1E5b", "lnav_l1l2": "L1L2"}
    user_type_name = None

    for pattern, type_name in user_type_name_def.items():
        if pattern in dset_name:
            return type_name

    if not user_type_name:
        log.fatal(f"User type name not found in {dset_name}.")


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

        log.info("Start: {} {}".format(program, " ".join(pandoc_args)))
        status = os.system(f"{program} {' '.join(pandoc_args)}")
        if status != 0:
            log.error("{} failed with error code {} ({})", program, status, " ".join([program] + pandoc_args))

        # TODO: pandoc subprocess call does not work. Why?
        # process = subprocess.run([program] + pandoc_args)
        # if process.returncode:
        #    log.error("{} failed with error code {} ({})", program, process.returncode, " ".join(process.args))


def _plot_bar_sisre_signal_combination_percentile(df, fid, figure_dir, threshold=False):
    """Generate bar plot with monthly SISRE 95% percentile for each GNSS signal combination

    Args:
       df (Dataframe):          Dataframe with time, satellite, system and GNSS signal combinations as columns
       fid (_io.TextIOWrapper): File object.
       figure_dir (PosixPath):  Figure directory.
       threshold (bool):        Plot threshold.
    """
    # df_monthly_percentile = df.set_index('time.gps').resample('D', how=lambda x: np.nanpercentile(x, q=95))
    df_monthly_percentile = df.dropna().set_index("time.gps").resample("M", how=lambda x: np.percentile(x, q=95))
    df_monthly_percentile.index = df_monthly_percentile.index.strftime("%b-%Y")
    df_monthly_percentile.transpose().plot(kind="bar")
    if threshold:
        plt.axhline(2, color="r")
        plot_filename = f"plot_bar_sisre_signal_combination_percentile_threshold.{FIGURE_FORMAT}"
    else:
        plot_filename = f"plot_bar_sisre_signal_combination_percentile.{FIGURE_FORMAT}"
    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISRE [m] (95th percentile)")
    plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3)
    plt.tight_layout()
    plt.savefig(figure_dir / plot_filename)
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly 95th percentile of global average SISRE for single- and dual-frequency users]({figure_dir}/{plot_filename})\n"
    )
    fid.write("\\newpage\n")


def _plot_bar_sisre_signal_combination_rms(df, fid, figure_dir):
    """Generate bar plot with monthly SISRE RMS for each GNSS signal combination

    Args:
       df (Dataframe):          Dataframe with time, satellite, system and GNSS signal combinations as columns
       fid (_io.TextIOWrapper): File object.
       figure_dir (PosixPath):  Figure directory.
    """
    df_monthly_rms = df.dropna().set_index("time.gps").resample("M", how=lambda x: np.sqrt(np.mean(np.square(x))))
    df_monthly_rms.index = df_monthly_rms.index.strftime("%b-%Y")
    df_monthly_rms.transpose().plot(kind="bar")
    plt.xlabel("Signal combination for single- and dual-frequency users")
    plt.xticks(rotation=0)
    plt.ylabel("SISRE RMS [m]")
    plt.legend(bbox_to_anchor=(1.04, 1), loc=2, borderaxespad=0., ncol=1)
    # plt.legend(bbox_to_anchor=(1, 1.15), loc=1, borderaxespad=0., ncol=3)
    plt.tight_layout()
    plt.savefig(figure_dir / f"plot_bar_sisre_signal_combination_rms.{FIGURE_FORMAT}")
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly RMS of global average SISRE for single- and dual-frequency users]({figure_dir}/plot_bar_sisre_signal_combination_rms.{FIGURE_FORMAT})\n"
    )
    fid.write("\\newpage\n")


def _plot_bar_sisre_satellite_percentile(df, fid, figure_dir, threshold=False):
    """Generate bar plot with monthly SISRE 95% percentile for each satellite

    Args:
       df (Dataframe):          Dataframe with time, satellite, system and GNSS signal combinations as columns
       fid (_io.TextIOWrapper): File object.
       figure_dir (PosixPath):  Figure directory.
    """

    user_types = ["E1", "E1E5a", "E1E5b"]  # TODO: better handling of user types !!!

    fig, axes = plt.subplots(len(user_types), 1, sharex=False, sharey=True)  # figsize=(6, 6));
    # plt.subplots_adjust(wspace=0.5, hspace=0.5);
    fig.set_figheight(20)  # inches

    if threshold:
        plot_filename = f"plot_bar_sisre_satellite_percentile_threshold.{FIGURE_FORMAT}"
    else:
        plot_filename = f"plot_bar_sisre_satellite_percentile.{FIGURE_FORMAT}"

    for ax, user in zip(axes, user_types):
        df_user = df.pivot(index="time.gps", columns="satellite", values=user)
        df_user_monthly_percentile = df_user.dropna().resample("M", how=lambda x: np.percentile(x, q=95))
        df_user_monthly_percentile.index = df_user_monthly_percentile.index.strftime("%b-%Y")
        df_user_monthly_percentile.transpose().plot(kind="bar", ax=ax, legend=False, title=user)
        ax.set_ylabel("SISRE [m] (95th percentile)")
        if threshold:
            plt.axhline(2, color="r")

    plt.legend(bbox_to_anchor=(0.8, -0.3), loc=1, borderaxespad=0., ncol=3)
    plt.tight_layout()
    plt.savefig(figure_dir / plot_filename)
    plt.clf()  # clear the current figure

    fid.write(
        f"![Monthly 95th percentile of global average SISRE of each satellite for single- and dual-frequency users]({figure_dir}/plot_bar_sisre_satellite_percentile.{FIGURE_FORMAT})\n"
    )
    fid.write("\\newpage\n")


def _header(fid):
    title = "Where"

    fid.write(
        """---
title: {title}
author: Where v{version} [{user}]
date: {nowdate:%Y-%m-%d}
---
""".format(
            # title=title["text"], version=where.__version__, user=config.analysis.user.str, nowdate=datetime.now()
            title=title,
            version=where.__version__,
            user=getpass.getuser(),
            nowdate=datetime.now(),  # TODO: Better solution?
        )
    )
