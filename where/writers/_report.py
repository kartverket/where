"""Class for reports

Description:
------------
TODO

"""
# Standard liberay imports
from datetime import date, datetime
import os
from shutil import copyfile
import subprocess
from typing import Union

# Where imports
import where
from where.lib import config
from where.lib import log
from where.lib import util


class Report:
    """Class for reports
    """

    def __init__(
        self, fid: "_io.TextIOWrapper", rundate: date, path: "pathlib.PosixPath", description: str = ""
    ) -> None:
        """Set up a new Report object

        Args:
            fid:         File object.
            rundate:     Run date.  
            path:        File path of markdown file.
            description: Short description of purpose of report (e.g. "RINEX navigation file analysis").
        """
        self.fid = fid
        self.rundate = rundate
        self.path = path
        self.description = description

    def markdown_to_pdf(self) -> None:
        """Convert markdown file to pdf format
        """
        # Close file object
        self.fid.close()

        if self.path.stat().st_size == 0:
            log.warn(f"Markdown file {self.path} is empty.")
            return 1

        pdf_path = str(self.path).replace(".md", ".pdf")
        program = "pandoc"

        # Convert markdown to pdf with pandoc
        pandoc_args = ["-f markdown", "-V classoption:twoside", "-N", "-o " + pdf_path, str(self.path)]
        log.info(f"Start: {program} {' '.join(pandoc_args)}")
        status = os.system(f"{program} {' '.join(pandoc_args)}")
        if status != 0:
            log.error(f"{program} failed with error code {status} ({' '.join([program] + pandoc_args)})")

        # TODO: pandoc subprocess call does not work. Why?
        # process = subprocess.run([program] + pandoc_args, shell=True)
        # if process.returncode:
        #    log.error(f"{program} failed with error code {process.returncode} ({' '.join(process.args)})")

    def write_config(self) -> None:
        """Print the configuration options for a given technique

        """
        # Print the individual configuration options
        self.fid.write(
            "#Configuration of {}\n__Located at {}__\n".format(self.description, ", ".join(config.tech.sources))
        )
        self.fid.write("```\n")
        self.fid.write(str(config.tech.as_str(key_width=25, width=70, only_used=True)))
        self.fid.write("\n```\n")
        self.fid.write("\\newpage\n\n")

    def write_dataframe_to_markdown(self, df: "Dataframe", format: str = "", statistic: bool = False) -> None:
        """Write Pandas DataFrame to Markdown table

        Args:
            df:         Pandas DataFrame
            format:     Define formatters for float columns (e.g. format = '6.3f')
            statistic:  Write statistical information rows at the end of the Markdown table
        """
        column_length = [len(c) for c in df.columns]

        # Check if dataframe is not empty
        if df.empty:
            log.warn("Nothing to write. Dataframe is empty.")
            return 0

        # Write header
        if list(df.index):  # Add DataFrame index to header
            num_space = len(max(df.index))
            head_line_1 = "\n| {} ".format(" " * num_space)
            head_line_2 = "|-{}-".format("-" * num_space)
        else:
            header_1 = ""

        self.fid.write(head_line_1 + "| {} |\n".format(" | ".join(list(df.columns))))
        self.fid.write(head_line_2 + "|-{}:|\n".format("-|-".join([n * "-" for n in column_length])))

        # Write data
        for index, row in df.iterrows():

            line = "| {idx:s} |".format(idx=index) if index else ""  # Add DataFrame index column

            for _, v in row.items():
                if isinstance(v, float):
                    line = line + " {:{fmt}} |".format(v, fmt=format)
                else:
                    line = line + " {} |".format(v)

            self.fid.write(line + "\n")

        # Write statistical information rows at the end of table
        if statistic:
            line_max = "| **max**  |"
            line_min = "| **min**  |"
            line_mean = "| **mean** |"
            for index, col in df.iteritems():
                line_max = line_max + " {:{fmt}} |".format(col.max(), fmt=format)
                line_min = line_min + " {:{fmt}} |".format(col.min(), fmt=format)
                line_mean = line_mean + " {:{fmt}} |".format(col.mean(), fmt=format)

            self.fid.write(line_max + "\n")
            self.fid.write(line_min + "\n")
            self.fid.write(line_mean + "\n")

        self.fid.write("\n")

    def add_figure(
        self, figure_path: Union["pathlib.PosixPath", str], caption: str = "", clearpage: bool = False
    ) -> None:
        """Add figure to report

        Args:
            figure_path:  Figure path.
            caption:      Caption of figure.
            clearpage:    Clear page.
        """
        self.fid.write(f"![{caption}]({figure_path})\n\n")

        if clearpage:
            self.fid.write("\n\\clearpage\n\n")

    def add_text(self, text: str, clearpage: bool = False) -> None:
        """Add text to report

        Args:
            text:       Text, which should be added to report
            clearpage:  Clear page.
        """
        self.fid.write(text + "\n\n")

        if clearpage:
            self.fid.write("\n\\clearpage\n\n")

    def title_page(self):
        """Write title page
        """

        title = f"{self.description} for day {self.rundate:%Y-%m-%d}"
        user_info = util.get_user_info()

        self.fid.write(
            """---
title: {title}
author: Where v{version} [{user}]
date: {nowdate:%Y-%m-%d}
---
    """.format(
                title=title,
                version=where.__version__,
                user=user_info.get("name", user_info["user"]),
                nowdate=datetime.now(),
            )
        )
        self.fid.write("\n\\newpage\n\n")
