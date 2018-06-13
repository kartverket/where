#!/usr/bin/env python3
"""Plot results of a Where analysis

Usage::

    {exe:there} [date] [--<pipeline>] [--session=<session>] [options]

The following commands are optional:

===================  ===========================================================
Command              Description
===================  ===========================================================
date                 The model run date in the format ``<year month day>``.
===================  ===========================================================

Furthermore, the following options are recognized:

===================  ===========================================================
Option               Description
===================  ===========================================================
{pipelines_doc:Plot results from}
--session=session    Plot results from the given session.

-S, --showconfig     Show configuration of There and exit.
--style=style        Set style of plotting. Choose between white, dark,
                     whitegrid, darkgrid, and ticks. Default is darkgrid.
--context=size       Set size of font and markers. Choose between paper,
                     notebook, talk, and poster. Default is notebook.
--colormap=cmap      Set colors used for plotting. Choose between all colormaps
                     recognized by matplotlib. Try --colormap=help for a list.
--debug, ...         Show additional debug information. Other flags such as
                     --all, --debug, --info, --warn, --error, --fatal, --none
                     are also allowed, and will show differing amounts of
                     information as the program runs.
--version            Show version information and exit.
-h, --help           Show this help message and exit.
===================  ===========================================================


Description:
------------

The script plots the results from a Where analysis for a given pipeline for the
given model run date. The analysis should already have been run with the main
where script.


Current Maintainers:
--------------------

{maintainers}

Version: {version}

"""

# Standard library imports
from datetime import datetime
import itertools
import os.path
import subprocess
import sys
import tkinter as tk
from tkinter import ttk
import webbrowser

# External library imports
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib
from mpl_toolkits.mplot3d import Axes3D  # noqa (Needed for projection='3d')
from matplotlib import cm
from matplotlib import gridspec
import numpy as np

# Where imports
from midgard.dev import console

import where
from where.lib import config
from where import data
from where.lib import exceptions
from where.lib import files
from where.lib import log
from where import pipelines
from where import setup
from where.lib.time import Time
from where.lib.unit import unit
from where.lib import util

# Optional import of editor and seaborn, add dummy methods in case they are not installed
from midgard.dev import optional

editor = optional.optional_import("editor", attrs=dict(edit=files.print_file))
sns = optional.optional_import("seaborn", attrs=dict(set=lambda **_: None))

PLOT_TYPES = dict()  # Dynamically set by plot_type()
matplotlib.rcParams["axes.formatter.useoffset"] = False


class CallWrapperReportingException(tk.CallWrapper):
    """Overrides the built-in CallWrapper, so that exceptions happening inside tkinter are reported."""

    def __call__(self, *args):
        """Apply first function SUBST to arguments, then FUNC."""
        try:
            if self.subst:
                args = self.subst(*args)
            return self.func(*args)
        except SystemExit:
            raise
        except:  # noqa
            import traceback
            from tkinter import messagebox

            err = traceback.format_exc()
            print(err)
            messagebox.showerror("Oops ...", err)


tk.CallWrapper = CallWrapperReportingException


@util.no_traceback
def main():
    """Parse command line arguments and run the program
    """
    # Use options to specify tech, include tech as profile
    # tech = None
    # for option, tech_name in pipelines.options().items():
    #     if util.check_options(option):
    #         tech = tech_name
    #         break
    try:
        tech = pipelines.get_from_options()
    except exceptions.UnknownPipelineError:
        tech = None

    session = util.read_option_value("--session", default="")

    # Add command line options to config
    user = config.program_vars(rundate=None, tech_name=None, session=None)["user"]
    profile = util.read_option_value("--profile")
    config.there.update_from_options(profile="__opts__", allow_new=True)
    config.there.profiles = [p for p in ("__opts__", profile, tech, user) if p]

    # Show the model configuration info
    if util.check_options("-S", "--showconfig"):
        print(config.there)
        return

    # Start logging
    log.init()
    if config.there.log_to_file.bool:
        log.file_init()
    else:
        log.turn_off_cache()

    # Initialize (poor man's optional commands)
    rundate = None
    util.parse_args(doc_module=__name__)
    if [a for a in sys.argv[1:] if not a.startswith("-")]:
        rundate = util.parse_args("date")

    # Use Seaborn to add styling
    sns.set(context=config.there.context.str, style=config.there.style.str)

    # Test if the colormap is recognized
    colormap = config.there.colormap.str
    if colormap:
        if colormap == "help":
            print("Possible values for colormap are:")
            print(console.fill(", ".join(sorted(cm.cmap_d)), initial_indent=" " * 4, subsequent_indent=" " * 4))
            sys.exit(0)
        else:
            cm.get_cmap(colormap)

    # Run program
    there = There(rundate, tech, session)
    there.mainloop()


def plot_type(func):
    """Decorator to register different plot types
    """
    name = func.__name__.replace("plot_", "").replace("_", " ")
    PLOT_TYPES[name] = func
    log.debug("Registered plot type '{}' handled by {}", name, func.__name__)

    return func


class There(tk.Tk):
    """A simple GUI wrapper around a Matplotlib chart for manipulating data

    The GUI makes it easy to interactively choose which data are plotted.
    """

    def __init__(self, rundate=None, tech=None, session=None, *args, **kwargs):
        """Set up the basic structure of the GUI
        """
        super().__init__(*args, **kwargs)
        super().wm_title(util.get_program_name().title())

        # Set up default values for variables
        self.vars = dict(
            config.program_vars(rundate, tech, session=session, location="work"), **config.date_vars(rundate)
        )
        if rundate is not None:
            self.vars.update({"date": f"{rundate:%Y%m%d}{session}{self.vars['id'].replace('_', '/', 1)}"})
        self.vars.update(config.there.initial_settings.as_dict())
        self.widget_updates = list()

        # Dataset controls
        dataset_line = ttk.Frame()
        self.add_dropdown(dataset_line, DD_Location)
        self.add_dropdown(dataset_line, DD_User)
        self.add_dropdown(dataset_line, DD_Date)
        self.add_dropdown(dataset_line, DD_Pipeline)
        self.add_dropdown(dataset_line, DD_Dataset)
        dataset_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Figure area
        self.figure = Plot(self, self.vars)
        self.add_widget_update(self.figure, self.figure.update_dataset)

        # Plot controls
        plot_line = ttk.Frame()
        self.add_dropdown(plot_line, DD_PlotType, self.figure)
        self.add_dropdown(plot_line, DD_XAxis, self.figure)
        self.add_dropdown(plot_line, DD_YAxis, self.figure)
        self.add_dropdown(plot_line, DD_Color, self.figure)
        self.add_dropdown(plot_line, DD_Size, self.figure)
        plot_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Filter controls
        filter_line = ttk.Frame()
        for filter in config.there.filters.list:
            self.add_dropdown(filter_line, filter_factory(filter), self.figure)
        filter_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Tool controls
        tool_line = ttk.Frame()
        self.add_checkbox(tool_line, "Scale to filter")
        self.add_checkbox(tool_line, "Force rerun")

        buttons = config.there.buttons.list
        for button in reversed(buttons):
            try:
                button_func = getattr(self, f"button_{button}")
            except AttributeError:
                available_buttons = ", ".join(b[7:] for b in dir(self) if b.startswith("button_"))
                raise ValueError(
                    f"Button '{button}' in config is not recognized. Available buttons are: {available_buttons}"
                ) from None
            self.add_button(tool_line, button.capitalize().replace("_", " "), button_func)
        tool_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Edit operations
        tool_line = ttk.Frame()
        double_click_opts = [o.capitalize().replace("_", " ") for o in config.there.double_click.list]
        self.add_radiobuttons(tool_line, "Double click", *double_click_opts)  # TODO: Remove event, ...

        # TODO: Generalize? Automatically make all ignore filter buttons and functions based on filters?
        self.add_button(tool_line, "Ignore baseline", self.figure.ignore_baseline)
        self.add_button(tool_line, "Ignore station", self.figure.ignore_station)
        self.add_button(tool_line, "Ignore source", self.figure.ignore_source)
        tool_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Update figure
        self.add_widget_update(self.figure, self.figure.update_plot)
        self.update_all()

    def add_button(self, frame, text, command):
        ttk.Button(frame, text=text, command=command).pack(side=tk.RIGHT, padx=1, pady=0)

    def add_checkbox(self, frame, text):
        CheckBox(frame, text=text, vars_=self.vars, figure=self.figure).pack(side=tk.LEFT, padx=5)

    def add_dropdown(self, frame, dd_cls, figure=None):
        tk.Label(frame, text=dd_cls.name.replace("_", " ").title() + ":").pack(side=tk.LEFT, padx=5, pady=5)
        dropdown = dd_cls(frame, self.vars, figure)
        dropdown.pack(side=tk.LEFT)
        self.add_widget_update(dropdown, dropdown.update_options)

    def add_radiobuttons(self, frame, group, *text_list):
        tk.Label(frame, text=group + ":").pack(side=tk.LEFT, padx=5)
        for text in text_list:
            Radiobutton(frame, text=text, group=group.replace(" ", "_").lower(), vars_=self.vars).pack(
                side=tk.LEFT, padx=5
            )

    def add_widget_update(self, widget, update_func):
        if self.widget_updates:
            last_widget = self.widget_updates[-1][0]
            last_widget.set_next_update(update_func)
        self.widget_updates.append((widget, update_func))

    def button_log(self):
        log_path = files.path("log", file_vars=self.vars)
        if not log_path.exists():
            log.warn("Log file {} does not exist", log_path)
            return

        log.print_file(log_path, config.there.print_log_level.str.upper())

    def button_show_map(self):
        map_path = files.path("output_web_map", file_vars=self.vars)
        if not map_path.exists():
            from where.writers import web_map

            web_map.web_map_writer(self.figure.dataset)
        webbrowser.open(map_path.as_uri())

    def button_config(self):
        cfg_path = files.path("config", file_vars=self.vars)
        if not cfg_path.exists():
            log.warn(f"Config file '{cfg_path}' does not exist")
            return

        editor.edit(cfg_path)
        setup.add_timestamp(self.vars["rundate"], self.vars["tech"], self.vars["session"], "last update")

    def button_rerun(self):
        """Rerun the current analysis

        Use subprocess to run the analysis, in order to have a separate process for config, log  etc
        """
        exe = where.__executable__
        cmd = "{exe} {rundate:%Y %-m %-d} --{tech} --session={dataset_name}".format(exe=exe, **self.vars).split()
        if self.vars["id"]:
            cmd.append(f"--id={self.vars['id'][1:]}")
        if self.vars["force_rerun"]:
            cmd.append("-F")
        log.info(f"Rerunning {self.vars['tech'].upper()} analysis: {' '.join(cmd)}")

        subprocess.check_call(cmd)
        self.update_all()

    def button_update(self):
        self.update_all()

    def button_remember(self):
        self.figure.remember_data()

    def button_forget(self):
        self.figure.forget_data()

    def update_all(self):
        self.widget_updates[0][1]()


class UpdateMixin():

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.next_update_func = None

    def set_next_update(self, next_update_func):
        self.next_update_func = next_update_func

    def update_next(self):
        if self.next_update_func is not None:
            self.next_update_func()


class Plot(FigureCanvasTkAgg, UpdateMixin):

    def __init__(self, master, vars_):
        figsize = (8, config.there.minimum_figure_height.float)
        self.figure = matplotlib.figure.Figure(figsize=figsize, dpi=100)
        self.graph = self.figure.add_subplot(1, 1, 1)
        super().__init__(self.figure, master)

        self.filters = list()
        self.vars = vars_
        self.dataset = None
        self.other_dataset = None
        self.next_update_func = None
        self.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        NavigationToolbar2TkAgg(self, master).update()
        self._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.picker_id = None
        self.cmap = config.there.colormap.str
        self.event_colors = dict()
        self.draw()

    def add_filter(self, filter_):
        self.filters.append(filter_.name)

    def update_dataset(self):
        if self.vars["rundate"] is None:
            location = config.files[f"directory_{self.vars['location']}"].directory.replaced.str.split("{")[0]
            log.error(f"No data found in {location}. Run '{where.__executable__}' to generate data.")
            return

        # Read dataset from disk
        self.dataset = data.Dataset(use_options=False, **self.vars)
        log.info("Updated: {}", self.dataset.repr)

        # Add event interval field
        events = self.dataset.get_events()
        if events and "time" in self.dataset.fields:
            obs_epochs = self.dataset.time.utc.jd
            event_epochs = np.array([e[0].utc.jd for e in self.dataset.get_events()])
            event, obs = np.meshgrid(event_epochs, obs_epochs)
            if "event_interval" not in self.dataset.fields:
                self.dataset.add_float("event_interval")
            self.dataset.event_interval[:] = np.sum(event - obs < 0, axis=1) + 0.5

        # Update the next widget
        self.update_next()

    def remember_data(self):
        self.other_dataset = self.dataset
        self.vars["x-axis_other"] = self.vars.get("x-axis_data")
        self.vars["y-axis_other"] = self.vars.get("y-axis_data")
        self.vars["color_other"] = self.vars.get("color_data")
        self.vars["size_other"] = self.vars.get("size_data")
        self.update_plot()

    def forget_data(self):
        self.other_dataset = None
        self.vars["x-axis_other"] = None
        self.vars["y-axis_other"] = None
        self.vars["color_other"] = None
        self.vars["size_other"] = None
        self.update_plot()

    def add_event_color(self, event_type):
        if event_type in self.event_colors:
            return

        event_types = list(self.event_colors.keys()) + [event_type]
        for event, value in zip(event_types, np.linspace(0, 1, len(event_types))):
            self.event_colors[event] = cm.tab20(value)

    @property
    def title(self):
        return "{} {}".format(self.vars["tech"].upper(), self.vars["date"])

    @property
    def xlabel(self):
        unit = self.dataset.unit(self.vars["x-axis_name"])
        unit_str = " [{}]".format(unit) if unit else ""
        return self.vars["x-axis_name"].replace("_", " ").replace(".", " - ").title() + unit_str

    @property
    def xlim(self):
        return self._calculate_range(self.vars["x-axis_data"], self.vars["x-axis_name"], "x-axis")

    @property
    def ylabel(self):
        unit = self.dataset.unit(self.vars["y-axis_name"])
        unit_str = " [{}]".format(unit) if unit else ""
        return self.vars["y-axis_name"].replace("_", " ").replace(".", " - ").title() + unit_str

    @property
    def ylim(self):
        return self._calculate_range(self.vars["y-axis_data"], self.vars["y-axis_name"], "x-axis")

    def update_plot(self):
        log.debug("Updating the {}-plot", self.vars["plot_type"])
        tooltip_fields = config.there.tooltip_fields.tuple
        idx, idx_other = self.do_filter()

        # Use the registered plotting functions to plot the correct plot
        self.figure.clear()
        PLOT_TYPES[self.vars["plot_type"]](self)

        # Some quick info, ignore if info can't be calculated
        info_data = self.vars["y-axis_data"]
        try:
            info_idx = idx & np.all(np.isfinite(info_data), axis=tuple(range(1, info_data.ndim)))
            log.info(
                "Num obs: {}  Mean: {}  RMS: {}",
                np.sum(info_idx),
                np.mean(info_data[info_idx]),
                np.sqrt(np.mean(np.square(info_data[info_idx]))),
            )
        except TypeError:
            pass

        if self.picker_id:
            self.figure.canvas.mpl_disconnect(self.picker_id)

        def on_pick(event):
            try:
                x_data = self.vars["x-axis_data"][idx]
                idx_finite_x = np.all(np.isfinite(x_data), axis=tuple(range(1, x_data.ndim)))
            except TypeError:
                idx_finite_x = np.ones(sum(idx), dtype=bool)
            try:
                y_data = self.vars["y-axis_data"][idx]
                idx_finite_y = np.all(np.isfinite(y_data), axis=tuple(range(1, y_data.ndim)))
            except TypeError:
                idx_finite_y = np.ones(sum(idx), dtype=bool)
            i_obs = np.where(idx)[0][idx_finite_x & idx_finite_y]

            pick_fields = [
                f
                for f in self.dataset.fields
                if f in tooltip_fields or any(f.startswith(p + "_") for p in tooltip_fields)
            ]
            fields = sorted(
                {
                    f
                    for f in pick_fields
                    + [
                        self.vars["x-axis_name"],
                        self.vars["y-axis_name"],
                        self.vars["color_name"],
                        self.vars["size_name"],
                    ]
                    if f in self.dataset.fields
                }
            )

            for ind in event.ind:
                texts = [
                    "{}: {}".format(f, self.dataset.plot_values(f)[i_obs[ind]]) for f in fields if f.startswith("time")
                ]
                texts += ["{}: {}".format(f, self.dataset[f][i_obs[ind]]) for f in fields if not f.startswith("time")]
                log.info("\n       ".join(texts))

        self.picker_id = self.figure.canvas.mpl_connect("pick_event", on_pick)
        self.figure.subplots_adjust(right=0.99, top=0.95)
        self.draw()

    def do_filter(self):
        filter_dict = {f: self.vars[f] for f in self.filters if self.vars[f] != "no filter"}
        idx_data = self.dataset.filter(**filter_dict)
        try:
            idx_other = self.other_dataset.filter(**filter_dict)
        except AttributeError:
            idx_other = idx_data
        return idx_data, idx_other

    def dbl_click_pick(self, event, mouse_event):
        if not mouse_event.dblclick:
            return False, dict(ind=list())

        dbl_click_func = getattr(self, f"dbl_click_{self.vars['double_click']}")
        return dbl_click_func(event, mouse_event)

    def dbl_click_do_nothing(self, _, __):
        log.debug("Doing nothing about double click")
        return False, dict()

    def dbl_click_add_clock_break(self, _, mouse_event):
        if self.vars["station"] == "no filter":
            log.error("Choose a station to add a clock break")
        else:
            time = Time(mouse_event.xdata, format="plot_date", scale="utc")

            # Check if there is a suspected clock break nearby
            all_data = self.vars["x-axis_data"]
            threshold = (max(all_data) - min(all_data)).total_seconds() * unit.seconds2days / 200
            for suspect_time, _, suspect_station in self.dataset.get_events("suspected_clock_break"):
                if suspect_station == self.vars["station"] and np.abs(time.utc.jd - suspect_time.utc.jd) < threshold:
                    log.info("Converting suspected clock break")
                    time = suspect_time
                    # TODO: dset.remove_event ...

            clock_break = f"{self.vars['station']} {time.datetime:{config.FMT_datetime}}"
            log.info(f"Adding clock break: '{clock_break}'")

            # Add event on dataset to visualize
            self.dataset.add_event(time, "unprocessed_clock_break", self.vars["station"])
            self.dataset.write()
            self.update_plot()

            # Add to config file
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                print(cfg.sources)
                current = cfg.vlbi_clock_correction.clock_breaks.as_list(", *")
                updated = ", ".join(sorted(current + [clock_break]))
                cfg.update("vlbi_clock_correction", "clock_breaks", updated, source=util.get_program_name())

        return False, dict()

    def ignore_baseline(self):
        if self.vars["baseline"] == "no filter":
            log.error("Choose a baseline in the filter menu to ignore it")
        else:
            log.info(f"Adding {self.vars['baseline']} to ignore_baseline")
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.vlbi_ignore_baseline.baselines.as_list(", *")
                updated = ", ".join(sorted(current + [self.vars["baseline"]]))
                cfg.update("vlbi_ignore_baseline", "baselines", updated, source=util.get_program_name())

    def ignore_station(self):
        if self.vars["station"] == "no filter":
            log.error("Choose a station in the filter menu to ignore it")
        else:
            log.info(f"Adding {self.vars['station']} to ignore_station")
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.ignore_station.stations.as_list(", *")
                updated = ", ".join(sorted(current + [self.vars["station"]]))
                cfg.update("ignore_station", "stations", updated, source=util.get_program_name())

    def ignore_source(self):
        if self.vars["source"] == "no filter":
            log.error("Choose a source in the filter menu to ignore it")
        else:
            log.info(f"Adding {self.vars['source']} to ignore_source")
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.vlbi_ignore_source.sources.as_list(", *")
                updated = ", ".join(sorted(current + [self.vars["source"]]))
                cfg.update("vlbi_ignore_source", "sources", updated, source=util.get_program_name())

    # Different types of plots
    @plot_type
    def plot_scatter(self, gs=gridspec.GridSpec(1, 1)[0]):
        idx, idx_other = self.do_filter()
        x_data, x_events = self._identify_events(self.vars["x-axis_data"])
        y_data, y_events = self._identify_events(self.vars["y-axis_data"])

        if x_data.ndim < 1 or y_data.ndim < 1:
            return idx, idx_other

        def event_pick(event, axis):
            all_data = self.vars[axis + "-axis_data"]
            try:
                threshold = (max(all_data) - min(all_data)).total_seconds() * unit.seconds2days / 100
            except AttributeError:
                threshold = (max(all_data) - min(all_data)) / 100
            e_time, e_type, e_description = event

            def on_pick(_, mouse_event):
                mouse_time = getattr(mouse_event, axis + "data")
                if mouse_time and abs(e_time.plot_date - mouse_time) < threshold:
                    log.info(
                        "Event: {} - {}\n       {}", e_time.datetime, e_type.replace("_", " ").title(), e_description
                    )
                    return True, dict(ind=list())
                else:
                    return False, dict(ind=list())

            return on_pick

        # Handle multi-dimensional data as separate plots
        ncols = 1 if x_data.ndim <= 1 else x_data.shape[1]
        nrows = 1 if y_data.ndim <= 1 else y_data.shape[1]
        sub_gs = gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=gs)

        for plot_num, (num_y, num_x) in enumerate(itertools.product(range(nrows), range(ncols))):
            ax = self.figure.add_subplot(sub_gs[plot_num])
            ax.clear()
            ax.scatter(0, 0, s=0, picker=self.dbl_click_pick, cmap=self.cmap)
            idx_x = slice(None) if x_data.ndim == 1 else (slice(None), num_x)
            idx_y = slice(None) if y_data.ndim == 1 else (slice(None), num_y)
            xlim = self.xlim if x_data.ndim == 1 else self.xlim[num_x]
            ylim = self.ylim if y_data.ndim == 1 else self.ylim[num_y]
            try:
                ax.scatter(
                    self.vars.get("x-axis_other")[idx_other][idx_x],
                    self.vars.get("y-axis_other")[idx_other][idx_y],
                    c=config.there.scatter.color_remember.str,
                    s=self.vars.get("size_other")[idx_other],
                    marker=config.there.scatter.marker_remember.str,
                    alpha=config.there.scatter.alpha_remember.float,
                )
            except (IndexError, TypeError, ValueError):
                pass
            ax.scatter(
                x_data[idx][idx_x],
                y_data[idx][idx_y],
                c=self._project_to_1d(self.vars["color_data"][idx]),
                s=self.vars["size_data"][idx],
                marker=config.there.scatter.marker.str,
                alpha=config.there.scatter.alpha.float,
                cmap=self.cmap,
                picker=True,
            )

            # Plot events
            for x in x_events:
                ax.plot(
                    (x[0].datetime, x[0].datetime),
                    self._pad_range(ylim),
                    ":",
                    color=self.event_colors[x[1]],
                    picker=event_pick(x, "x"),
                )
            for y in y_events:
                ax.plot(
                    self._pad_range(xlim),
                    (y[0].datetime, y[0].datetime),
                    ":",
                    color=self.event_colors[y[1]],
                    picker=event_pick(y, "y"),
                )

            # Label subplot
            ax.set_xlim(self._pad_range(xlim))
            ax.set_ylim(self._pad_range(ylim))
            if num_x > 0:
                ax.tick_params(left=False, labelleft=False)
            if num_y < nrows - 1:
                ax.tick_params(bottom=False, labelbottom=False)
            else:
                for label in ax.get_xticklabels():
                    label.set_rotation(12)

        # Label figure
        self.figure.suptitle(self.title)
        self.figure.text(0.5, 0.01, self.xlabel, ha="center")
        self.figure.text(0.01, 0.5, self.ylabel, va="center", rotation="vertical")

        # Return indices to surrounding plot functions
        return idx, idx_other

    @plot_type
    def plot_scatter_w_hist(self):
        gs = gridspec.GridSpec(2, 2, height_ratios=[1, 5], width_ratios=[5, 1], hspace=.02, wspace=.01)
        idx, idx_other = self.plot_scatter(gs=gs[2])

        # Horisontal histogram
        histogram = self.figure.add_subplot(gs[0])
        hist_data, hist_range = self._convert_datetime_to_sec(self.vars["x-axis_data"])
        hist_range = self._pad_range(hist_range)

        histogram.hist(hist_data[idx], bins=99, range=hist_range, alpha=0.75)
        histogram.set_xlim(hist_range)
        histogram.tick_params(bottom=False, labelbottom=False)
        histogram.set_ylabel("Count")

        # Vertical histogram
        histogram = self.figure.add_subplot(gs[3])
        hist_data, hist_range = self._convert_datetime_to_sec(self.vars["y-axis_data"])
        hist_range = self._pad_range(hist_range)

        histogram.hist(hist_data[idx], bins=49, range=hist_range, alpha=0.75, orientation="horizontal")
        histogram.set_ylim(hist_range)
        histogram.tick_params(left=False, labelleft=False)
        histogram.set_xlabel("Count")

    @plot_type
    def plot_polar(self):
        ax = self.figure.add_axes([0.1, 0.05, 0.8, 0.8], polar=True)
        ax.clear()

        idx, idx_other = self.do_filter()
        x_data, x_lim = self._convert_datetime_to_sec(self._project_to_1d(self.vars["x-axis_data"]))
        y_data, y_lim = self._convert_datetime_to_sec(self._project_to_1d(self.vars["y-axis_data"]))

        # Convert radians to degrees for certain fields
        y_name = self.vars["y-axis_name"]
        if y_name.endswith(".zenith_distance") or y_name.endswith(".elevation") or y_name.endswith(".azimuth"):
            # Todo: Use unit instead to convert radians to degrees
            y_data = np.degrees(y_data)
            y_lim = (np.min((0, np.degrees(np.min(y_lim)))), np.max((90, np.degrees(np.max(y_lim)))))

        ax_scatter = ax.scatter(
            x_data[idx],
            y_data[idx],
            c=self.vars["color_data"][idx],
            s=self.vars["size_data"][idx],
            marker="x",
            alpha=0.7,
            cmap=self.cmap,
            picker=True,
        )

        ax.set_ylim(self._pad_range(y_lim))
        ax.set_theta_zero_location("N")  # sets 0(deg) to North
        ax.set_theta_direction(-1)  # sets plot clockwise
        if self.vars["color_name"]:
            self.figure.colorbar(ax_scatter)

        ax.set_title("{} vs {}".format(self.xlabel, self.ylabel))
        self.figure.suptitle(self.title)

    @plot_type
    def plot_3d(self):
        idx, idx_other = self.do_filter()

        # Pick out data: Prefer y-axis, but alternatively plot x-axis if it is 3d
        if self.vars["y-axis_3d"]:
            plot_data = self.vars["y-axis_data"]
            x_name, y_name, z_name = ["{}: {}".format(xyz, self.ylabel) for xyz in "XYZ"]
        elif self.vars["x-axis_3d"]:
            plot_data = self.vars["x-axis_data"]
            x_name, y_name, z_name = ["{}: {}".format(xyz, self.xlabel) for xyz in "XYZ"]
        else:
            data_x, _ = self._convert_datetime_to_sec(self.vars["x-axis_data"])
            data_y, _ = self._convert_datetime_to_sec(self.vars["y-axis_data"])
            plot_data = np.vstack((data_x, data_y, np.zeros(data_x.shape))).T
            x_name, y_name, z_name = self.xlabel, self.ylabel, ""

        # Do the plotting
        ax = self.figure.add_subplot(1, 1, 1, projection="3d")
        ax.clear()
        ax.scatter(
            plot_data[idx, 0],
            plot_data[idx, 1],
            plot_data[idx, 2],
            c=self.vars["color_data"][idx],
            s=self.vars["size_data"][idx],
            marker="x",
            cmap=self.cmap,
            picker=True,
        )

        ax.set_title(self.title)
        ax.set_xlabel(x_name)
        ax.set_ylabel(y_name)
        ax.set_zlabel(z_name)

    def _convert_datetime_to_sec(self, plot_data):
        if np.issubdtype(plot_data.dtype, np.datetime64) or plot_data.dtype == np.object:
            time_data = np.array([t.timestamp() for t in plot_data])
            plot_data = time_data - time_data[0]
        plot_range = self._calculate_range(plot_data.flatten(), "field", "plot data")

        return plot_data, plot_range

    def _identify_events(self, plot_data):
        events = self.dataset.get_events()
        for _, event_type, _ in events:
            self.add_event_color(event_type)

        if events and (np.issubdtype(plot_data.dtype, np.datetime64) or plot_data.dtype == np.object):
            return plot_data, events
        else:
            return plot_data, list()

    @staticmethod
    def _project_to_1d(plot_data):
        return plot_data.mean(axis=tuple(range(1, plot_data.ndim))) if plot_data.ndim > 1 else plot_data

    def _calculate_range(self, values, name, description):
        """Calculate range of given data
        """
        if self.vars["scale_to_filter"]:
            idx_filter, _ = self.do_filter()
        else:
            idx_filter = np.ones(len(values), dtype=bool)

        values = values[idx_filter]
        try:
            idx_finite = np.all(np.isfinite(values), axis=tuple(range(1, values.ndim)))
        except TypeError:  # Datetimes throw a TypeError on isnan/isfinite
            idx_finite = np.ones(len(values), dtype=bool)

        if not np.any(idx_finite):
            log.warn("No data for '{}' ({})", name, description)
            return 0, 1

        if values.ndim == 1:
            return np.min(values[idx_finite]), np.max(values[idx_finite])
        else:
            return tuple(zip(np.min(values[idx_finite], axis=0), np.max(values[idx_finite], axis=0)))

    @staticmethod
    def _pad_range(val_range, factor=0.01):
        """Pad value range at both sides
        """
        rmin, rmax = val_range
        try:
            delta = np.fmax((rmax - rmin) * factor, 1e-6)
        except TypeError:  # np.fmax fails on timedeltas
            delta = (rmax - rmin) * factor

        return rmin - delta, rmax + delta


class CheckBox(ttk.Checkbutton):

    def __init__(self, master, text, vars_, figure, **kwargs):
        self.choice = tk.IntVar()
        super().__init__(master, text=text, command=self.update_vars, variable=self.choice, **kwargs)
        self.name = text.replace(" ", "_").lower()
        self.vars = vars_
        self.figure = figure
        self.vars[self.name] = False

    def update_vars(self):
        self.vars[self.name] = bool(self.choice.get())
        self.figure.update_plot()
        log.debug("Setting {} = {}", self.name, self.vars[self.name])


class Radiobutton(ttk.Radiobutton):

    groupvars = dict()

    def __init__(self, master, text, group, vars_, **kwargs):
        self.name = text.replace(" ", "_").lower()
        if group not in self.groupvars:
            self.groupvars[group] = tk.StringVar()
            self.groupvars[group].set(self.name)
        self.choice = self.groupvars[group]
        super().__init__(master, text=text, command=self.update_vars, variable=self.choice, value=self.name, **kwargs)
        self.group = group
        self.vars = vars_
        self.update_vars()

    def update_vars(self):
        self.vars[self.group] = self.choice.get()
        log.debug("Setting {} = {}", self.group, self.vars[self.group])


class DropdownPicker(ttk.Combobox, UpdateMixin):

    name = "no_name"
    width = 8

    def __init__(self, master, vars_, figure, **kwargs):
        self.choice = tk.StringVar()
        super().__init__(master, width=self.width, textvariable=self.choice, **kwargs)

        self.vars = vars_
        self.figure = figure
        self.file_vars = config.files.vars  # TODO: This should not depend on config.files.vars
        self.next_update_func = None
        self.choice.set(self.vars.get(self.name, ""))
        self.choice.trace("w", self.choose)

    def _glob_dset_paths(self, file_vars):
        # TODO: How to deal with archive vs work
        # print(self.vars)
        # print(file_vars)
        ...

    def choose(self, *_):
        vars_ = self.parse_vars()
        self.vars.update(vars_)
        for name, value in vars_.items():
            log.debug("Setting {} = {}", name, value)
        self.update_next()

    def parse_vars(self):
        return {self.name: self.choice.get()}

    def update_options(self):
        log.debug("Updating {}", self.name)
        self["values"] = self.read_options()
        if self["values"] and self.choice.get() not in self["values"]:
            self.choice.set(self["values"][0])
        elif not self["values"]:
            self.choice.set("")
        else:
            self.choice.set(self.choice.get())

    def read_options(self):
        return list()


class DD_Location(DropdownPicker):

    name = "location"
    width = 8

    def read_options(self):
        """Read possible locations from the file config
        """
        return sorted([d[10:] for d in config.files.sections if d.startswith("directory_")])


class DD_User(DropdownPicker):

    name = "user"
    width = 8

    def read_options(self):
        """Read users from the file directories
        """
        users = files.glob_variable("dataset_hdf5", "user", r"[a-z]+")
        return sorted(users)


class DD_Date(DropdownPicker):

    name = "date"
    width = 30

    def read_options(self):
        """Read dates from filenames

        TODO: The persistence of config.files.vars causes problems. We are not able to glob for paths since date is
              already set in the vars.
        """
        file_vars = dict(user=self.vars["user"])
        paths = files.glob_paths("dataset_hdf5", file_vars=file_vars)
        dirs = [os.path.basename(os.path.dirname(p)) for p in paths]
        dates = set()
        for dirname in dirs:
            dates.add(dirname.replace("_", "/", 1))
        return sorted(dates)

    def parse_vars(self):
        """Construct date variables and split out timestamp
        """
        parts = (self.choice.get() + "/").split("/")
        rundate = datetime.strptime(parts[0][:8], "%Y%m%d").date() if parts[0] else None
        file_vars = config.date_vars(rundate) if rundate else dict()
        file_vars["rundate"] = rundate
        file_vars["date"] = parts[0][:8]
        file_vars["session"] = parts[0][8:]
        file_vars["timestamp"] = parts[1].split("_")[0]
        try:
            # Test if timestamp is valid
            datetime.strptime(file_vars["timestamp"], config.FMT_dt_file)
            parts[1] = "_".join(parts[1].split("_")[1:])
        except ValueError:
            file_vars["timestamp"] = ""
        file_vars["id"] = "_" + parts[1] if parts[1] else ""
        return file_vars


class DD_Pipeline(DropdownPicker):

    name = "pipeline"
    width = 14

    def read_options(self):
        """Read pipeline and stage from filenames
        """
        file_vars = {k: v for k, v in self.vars.items() if k not in ("tech", "stage")}
        paths = files.glob_paths("dataset_hdf5", file_vars=file_vars)
        tech_stages = [(s[0], s[2]) for s in [os.path.basename(p).split("-") for p in sorted(paths)]]
        return [f"{t}/{s}" for t, s in sorted(tech_stages, key=self._sorter)]

    @staticmethod
    def _sorter(tech_stage):
        """Return a sort key for the given tech and stage

        Sorts first on tech then on stage such that later stages are sorted first. Unknown stages are sorted
        alphabetically at the end. Stages can be prioritized by adding them to stages in there.conf.

        Args:
            tech_stage:  Tuple with tech and stage.

        Returns:
            Tuple that can be sorted on.
        """
        tech, stage = tech_stage
        stages = config.there.stages.list + pipelines.stages(tech)[::-1]  # Reversed to sort later stages first
        try:
            stage_id = stages.index(stage)
        except ValueError:
            stage_id = len(stages)  # Unknown stages are sorted last
        return (tech, stage_id, stage)

    def parse_vars(self):
        if not self.choice.get():
            return dict()

        tech, _, stage = self.choice.get().partition("/")
        return dict(tech=tech, stage=stage)


class DD_Dataset(DropdownPicker):

    name = "dataset"
    width = 14

    def read_options(self):

        """Read dataset name and id
        """
        return sorted(data.list_datasets(use_options=False, **self.vars), reverse=True)

    def parse_vars(self):
        dataset_name, _, dataset_id = self.choice.get().partition("/")
        try:
            return dict(dataset_name=dataset_name, dataset_id=int(dataset_id))
        except ValueError:
            return dict(dataset_name=dataset_name, dataset_id=0)


class DD_PlotType(DropdownPicker):

    name = "plot_type"
    width = 10

    def read_options(self):
        return sorted(PLOT_TYPES)


class DD_Field(DropdownPicker):

    width = 20

    def read_options(self):
        return sorted(self.figure.dataset.fields, key=lambda f: ("a" if f.startswith("time") else "z") + f)

    def parse_vars(self):
        field = self.choice.get()
        if not field or field == "none":
            values = np.ones(self.figure.dataset.num_obs)
            return {self.name + "_data": values, self.name + "_name": ""}

        values = self.figure.dataset.plot_values(field)
        do_3d_plot = values.ndim == 2 and values.shape[1] == 3

        # Handle pairs of 3d-values as differences
        # if values.ndim == 2 and values.shape[1] == 6:
        #    values = np.diff(values.reshape(len(values), -1, 3), axis=1)[:, 0, :]
        #    do_3d_plot = True

        return {self.name + "_data": values, self.name + "_name": field, self.name + "_3d": do_3d_plot}


class DD_XAxis(DD_Field):

    name = "x-axis"


class DD_YAxis(DD_Field):

    name = "y-axis"


class DD_Color(DD_Field):

    name = "color"

    def read_options(self):
        return ["none"] + super().read_options()

    def parse_vars(self):
        file_vars = super().parse_vars()
        try:
            file_vars[self.name + "_data"] = np.vectorize(lambda d: d.timestamp())(file_vars[self.name + "_data"])
        except AttributeError:
            pass

        return file_vars


class DD_Size(DD_Field):

    name = "size"

    def read_options(self):
        return ["none"] + super().read_options()

    def parse_vars(self):
        fields = super().parse_vars()
        values = fields[self.name + "_data"]
        try:
            values = np.vectorize(lambda d: d.timestamp())(values)
        except AttributeError:
            pass

        size = normalize(np.abs(values), self.figure)  # TODO: Can we normalize while plotting to update properly?
        return dict(size_data=size ** 2 * 65 + 10, size_name=fields[self.name + "_name"])


def filter_factory(filter_name):

    class DD_Filter(DropdownPicker):

        width = 20
        name = filter_name

        def __init__(self, master, file_vars, figure, **kwargs):
            super().__init__(master, file_vars, figure, **kwargs)
            figure.add_filter(self)

        def read_options(self):
            try:
                return ["no filter"] + self.figure.dataset.unique(self.name)
            except KeyError:
                return ["no filter"]

    return DD_Filter


def normalize(values, figure):
    """Normalize the values of values so they are in the [0, 1]-interval
    """
    if figure.vars["scale_to_filter"]:
        idx_filter, _ = figure.do_filter()
    else:
        idx_filter = np.ones(len(values), dtype=bool)

    min_ = np.min(values[idx_filter])
    max_ = np.max(values[idx_filter])
    rng = max_ - min_

    if rng == 0:
        return 0.5 * np.ones(values.shape)

    normalized_vals = (values - min_) / rng
    normalized_vals[~idx_filter] = 0.5

    return normalized_vals


# Run main function only when running as script
if __name__ == "__main__":
    main()
