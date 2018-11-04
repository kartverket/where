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
import calendar
from collections import namedtuple
from contextlib import contextmanager
from datetime import date, datetime, timedelta
import itertools
import subprocess
import sys
import tkinter as tk
from tkinter import ttk
from tkinter import scrolledtext
from tkinter import font as tk_font
import webbrowser

# External library imports
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
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
TABS = dict()  # Dynamically set by register_tab()
DROPDOWNS = dict()  # Dynamically set by register_dropdown()
KEY_EVENTS = dict()  # Read from config by set_up_key_events()
UPDATE = True  # Global flag

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
    # Use options to specify pipeline and session, include pipeline as profile
    try:
        pipeline = pipelines.get_from_options()
    except exceptions.UnknownPipelineError:
        pipeline = None
    session = util.read_option_value("--session", default="")  # Not using pipelines.get_session to avoid validation

    # Add command line options to config
    user = config.program_vars(rundate=None, tech_name=None, session=None)["user"]
    profile = util.read_option_value("--profile")
    config.there.update_from_options(profile="__opts__", allow_new=True)
    config.there.profiles = [p for p in ("__opts__", profile, pipeline, user) if p]

    # Show the model configuration info
    if util.check_options("-S", "--showconfig"):
        print(config.there)
        raise SystemExit

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
            raise SystemExit
        else:
            cm.get_cmap(colormap)

    # Run program
    there = There(rundate, pipeline, session)
    there.mainloop()


def plot_type(func):
    """Decorator to register different plot types
    """
    name = func.__name__.replace("plot_", "").replace("_", " ")
    PLOT_TYPES[name] = func
    log.debug("Registered plot type '{}' handled by {}", name, func.__name__)

    return func


def register_tab(cls):
    """Decorator to register tabs"""
    TABS[cls.name] = cls
    log.debug(f"Registered tab '{cls.name}' handled by {cls.__name__}")


def register_dropdown(cls):
    """Decorator to register dropdowns"""
    DROPDOWNS[cls.name] = cls
    log.debug(f"Registered dropdown '{cls.name}' handled by {cls.__name__}")


class There(tk.Tk):
    """A simple GUI wrapper around a Matplotlib chart for manipulating data

    The GUI makes it easy to interactively choose which data are plotted.
    """

    def __init__(self, rundate=None, pipeline=None, session=None, *args, **kwargs):
        """Set up the basic structure of the GUI
        """
        super().__init__(*args, **kwargs)
        super().wm_title(util.get_program_name().title())

        self._there_icon = tk.PhotoImage(file=files.path("there_icon"))  # Keep reference to icon
        self.iconphoto(True, self._there_icon)

        # Set up default values for variables
        self.vars = dict(
            common=dict(
                config.program_vars(rundate, pipeline, session=session), **config.date_vars(rundate), rundate=rundate
            )
        )
        if rundate is not None:
            self.vars["common"].update({"date_session": f"{rundate:%Y%m%d}/{session}"})
        self.vars["common"]["pipeline"] = self.vars["common"]["tech"]  # TODO: Remove when tech is changed to pipeline
        self.vars["common"]["do_update"] = True
        self.vars["common"]["id"] = self.vars["common"]["id"].lstrip("_")

        # Add tabs
        self.status = Status(self)
        self.nb = ttk.Notebook(self)
        self.nb.enable_traversal()
        self.tabs = dict()

        tabs = config.there.tabs.list
        for tab in tabs:
            # Tab variables
            self.vars[tab] = self.vars["common"].copy()
            if f"{tab}_initial" in config.there.section_names:
                self.vars[tab].update(config.there[f"{tab}_initial"].as_dict())

            # Add tab widget
            self.tabs[tab] = TABS[tab](self.nb, self.vars[tab])
            self.tabs[tab].add_to(self.nb)
        self.nb.pack(expand=True, fill="both")

        # Only display status line if set in config
        if config.there.status.display.bool:
            self.status.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Select tab
        if pipeline and rundate:
            self.select_tab("session")
        elif pipeline:
            self.select_tab("timeseries")
        else:
            self.select_tab("run_analysis")

        # Test key events
        self.set_up_key_events()

    def set_up_key_events(self):
        # Method events are first tried on the current tab, then on There in general
        if "key_method" in config.there.section_names:
            for key, entry in config.there.key_method.items():
                KEY_EVENTS[key] = (entry.str.capitalize().replace("_", " "), [["current_tab", entry.str], [entry.str]])

        # Config events add to the current tech config using update_config()
        if "key_config" in config.there.section_names:
            for key, entry in config.there.key_config.items():
                KEY_EVENTS[key] = (
                    f"Update config: {entry.str}",
                    [["current_tab", "update_config"]],
                    entry.as_list(", *"),
                )

        # Select events select from dropdowns using select_dropdown()
        if "key_select" in config.there.section_names:
            for key, entry in config.there.key_select.items():
                KEY_EVENTS[key] = (f"Select {entry.str}", [["current_tab", "select_dropdown"]], entry.as_list(", *"))

        self.bind("<Key>", self.key_pressed)

    def key_pressed(self, event):
        key = event.char
        if key not in KEY_EVENTS:
            return

        event, attrs_list, *args = KEY_EVENTS[key]
        for attrs in attrs_list:
            obj = self
            try:
                for attr in attrs:
                    obj = getattr(obj, attr)
                obj(*args)
                break
            except AttributeError:
                continue
        else:
            log.warn(f"Key {key!r} is not valid on the {self.current_tab.name.title().replace('_', ' ')!r} tab")

    @property
    def current_tab(self):
        tab = self.nb.select().rpartition(".")[-1]
        return self.tabs[tab]

    def select_tab(self, tab):
        if tab in self.tabs:
            return self.tabs[tab].select()
        log.warn(f"Could not select unknown tab {tab!r}")

    def show_keyboard_shortcuts(self):
        key_str = ["Current Keyboard Shortcuts:"]
        key_str += [f"{k}: {v[0]}" for k, v in sorted(KEY_EVENTS.items())]
        log.out("\n  ".join(key_str))
        self.current_tab.status("\n".join(key_str))


class Status(ttk.Frame):
    def __init__(self, master):
        super().__init__(master)

        cfg = config.there.status

        # Set up fonts
        font_family, _, _ = font = (cfg.font_family.str, cfg.font_size.int, cfg.font_style.str)
        if font_family:
            available_fonts = tk_font.families()
            if font_family not in available_fonts:
                log.warn(f"Font family {font_family!r} is not available. Use one of: {', '.join(available_fonts)}")
        else:
            font_family = tk_font.nametofont("TkTextFont").actual()["family"]
            font = (font_family,) + font[1:]

        # Lay out widgets
        self._previous_icon = tk.PhotoImage(file=files.path("there_button_icon", file_vars=dict(name="previous")))
        self._previous = ttk.Button(self, image=self._previous_icon, command=lambda: self.write("Previous"))
        # self._previous.pack(side=tk.LEFT, fill=tk.Y)
        self._next_icon = tk.PhotoImage(file=files.path("there_button_icon", file_vars=dict(name="next")))
        self._next = ttk.Button(self, image=self._next_icon, command=lambda: self.write("Next"))
        # self._next.pack(side=tk.RIGHT, fill=tk.Y)
        self._text = scrolledtext.ScrolledText(self, height=cfg.height.int, wrap=tk.WORD, font=font)
        self._text.pack(fill=tk.X, expand=False)

        # Initialize widgets
        startup_text = f"Start {util.get_program_info()} at {datetime.now():{config.FMT_datetime}}"
        self._text.insert(1.0, startup_text)
        self._text.config(state=tk.DISABLED)

        self._previous.config(state=tk.DISABLED)
        self._next.config(state=tk.DISABLED)

    def clear(self):
        self._text.config(state=tk.NORMAL)
        self._text.delete(1.0, tk.END)
        self._text.config(state=tk.DISABLED)

    def write(self, text):
        self._text.config(state=tk.NORMAL)
        self._text.insert(tk.END, f"\n{text}")  # New-line in front to avoid blank line at end
        self._text.yview(tk.END)  # Move scrollbar to bottom
        self._text.config(state=tk.DISABLED)


class Tab(ttk.Frame):

    name = "updated by subclasses"

    def __init__(self, master, vars):
        super().__init__(master, name=self.name)
        self.master = master
        self.vars = vars
        self.widgets = dict()
        self.widget_updates = list()
        self._icon = None
        self.init_tab()

    def init_tab(self):
        """Initialize the contents of the tab"""
        pass

    @property
    def icon(self):
        if self._icon is None:
            self._icon = tk.PhotoImage(file=files.path("there_tab_icon", file_vars=dict(name=self.name)))
        return self._icon

    def add_to(self, master):
        """Add tab to master widget"""
        master.add(
            self, text=f"  {self.name.title().replace('_', ' ')}  ", underline=2, image=self.icon, compound="left"
        )

    def select(self):
        self.master.select(f"{self.master}.{self.name}")

    def status(self, text):
        self.master.master.status.write(text)

    def add_button(self, frame, text, command):
        ttk.Button(frame, text=text, command=command).pack(side=tk.RIGHT, padx=1, pady=0)

    def add_buttons(self, frame, instance, *button_list):
        for button in reversed(button_list):  # Reversed since we are packing to the right
            try:
                button_func = getattr(instance, f"button_{button}")
            except AttributeError:
                buttons = ", ".join(b[7:] for b in dir(instance) if b.startswith("button_"))
                raise ValueError(f"Button '{button}' in config is not recognized. Use one of: {buttons}") from None
            self.add_button(frame, button.title().replace("_", " "), button_func)

    def add_checkboxes(self, frame, instance, *name_list):
        for name in name_list:
            CheckBox(frame, name=name, vars_=self.vars, figure=instance).pack(side=tk.LEFT, padx=5)

    def add_dropdown(self, frame, name, figure=None):
        try:
            dropdown = DROPDOWNS[name](frame, self.vars, figure)
        except KeyError:
            dropdowns = ", ".join(DROPDOWNS.keys())
            raise ValueError(f"Dropdown '{name}' in config is not recognized. Use one of: {dropdowns}") from None
        self.widgets[name] = dropdown
        tk.Label(frame, text=f"{dropdown.label}:").pack(side=tk.LEFT, padx=5, pady=5)
        dropdown.pack(side=tk.LEFT)
        self.add_widget_update(dropdown, dropdown.update_options)

    def add_entry(self, frame, name):
        tk.Label(frame, text=f"{name.title()}:").pack(side=tk.LEFT, padx=5, pady=5)
        entry = tk.Entry(frame, width=20)
        entry.pack(side=tk.LEFT)

    def add_widget_update(self, widget, update_func):
        if self.widget_updates:
            last_widget = self.widget_updates[-1][0]
            last_widget.set_next_update(update_func)
        self.widget_updates.append((widget, update_func))

    def add_radiobuttons(self, frame, group, *name_list):
        self.widgets[group] = dict()
        group_name = group.replace("_", " ").title()
        tk.Label(frame, text=f"{group_name}:").pack(side=tk.LEFT, padx=5)
        for name in name_list:
            self.widgets[group][name] = Radiobutton(frame, name=name, group=group, vars_=self.vars)
            self.widgets[group][name].pack(side=tk.LEFT, padx=5)

    def update_all(self):
        self.widget_updates[0][1]()

    def remember_data(self):
        return self.figure.remember_data()

    def forget_data(self):
        return self.figure.forget_data()

    def show_vars(self):
        from pprint import pprint

        pprint(self.vars)

    update_figure = update_all
    button_update = update_figure
    button_remember = remember_data
    button_forget = forget_data
    button_show_vars = show_vars

    def button_run_analysis(self):
        """Run the current analysis

        Use subprocess to run the analysis, in order to have a separate process for config, log  etc
        """
        exe = where.__executable__
        cmd = "{exe} {rundate:%Y %-m %-d} --{pipeline} --session={dataset_name}".format(exe=exe, **self.vars).split()
        if self.vars["id"]:
            cmd.append(f"--id={self.vars['id'][1:]}")
        if self.vars.get("force_rerun"):
            cmd.append("-F")
        log.info(f"Running {self.vars['pipeline'].upper()} analysis: {' '.join(cmd)}")

        subprocess.check_call(cmd)
        self.update_all()

    def button_rerun(self):
        """Rerun the current analysis

        Use subprocess to run the analysis, in order to have a separate process for config, log  etc
        """
        exe = where.__executable__
        cmd = "{exe} {rundate:%Y %-m %-d} --{pipeline} --session={dataset_name}".format(exe=exe, **self.vars).split()
        if self.vars["id"]:
            cmd.append(f"--id={self.vars['id'][1:]}")
        if self.vars.get("force_rerun", False):
            cmd.append("-F")
        log.info(f"Rerunning {self.vars['pipeline'].upper()} analysis: {' '.join(cmd)}")

        subprocess.check_call(cmd)
        self.update_all()

    def select_dropdown(self, selections):
        for selection in selections:
            key, _, value = selection.partition("=")
            if key not in self.widgets:
                continue
            if value not in self.widgets[key]["values"]:
                continue
            self.widgets[key].choice.set(value)


class TabFigure(Tab):
    """General figure tab"""

    def init_tab(self):
        cfg = config.there[self.name]

        # Dataset controls
        dataset_line = ttk.Frame(self)
        for dropdown in cfg.dataset_dropdowns.list:
            self.add_dropdown(dataset_line, dropdown)
        dataset_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Figure area
        self.figure = Plot(self, self.vars)
        self.add_widget_update(self.figure, self.figure.update_dataset)

        # Plot controls
        plot_line = ttk.Frame(self)
        self.add_dropdown(plot_line, "plot_type", self.figure)  # TODO: Config?
        self.add_dropdown(plot_line, "x_axis", self.figure)
        self.add_dropdown(plot_line, "y_axis", self.figure)
        self.add_dropdown(plot_line, "color", self.figure)
        self.add_dropdown(plot_line, "size", self.figure)
        plot_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Filter controls
        filter_line = ttk.Frame(self)
        for filter in cfg.filters.list:
            filter_cls = filter_factory(filter)
            self.add_dropdown(filter_line, filter_cls.name, self.figure)
        filter_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Tool controls
        button_line = ttk.Frame(self)
        self.add_checkboxes(button_line, self.figure, *cfg.checkboxes.list)
        self.add_buttons(button_line, self, *cfg.buttons.list)
        button_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Edit operations
        tool_line = ttk.Frame(self)
        self.add_radiobuttons(tool_line, "double_click", *cfg.double_clicks.tuple)
        self.add_buttons(tool_line, self.figure, *cfg.figure_buttons.list)
        tool_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Update figure
        self.add_widget_update(self.figure, self.figure.update_plot)
        self.update_all()

    def dataset_info(self):
        dset = self.figure.dataset
        if dset is None:
            return

        update_txt = [dset.repr]
        timestamps = config.timestamps(**self.vars)
        if timestamps:
            update_txt += [f"{k.title()} {v}" for k, v in timestamps.items()]
        log.out("\n  ".join(update_txt))
        self.status(", ".join(update_txt))

    def next_double_click_option(self):
        self.next_radiobutton_option("double_click")

    def next_radiobutton_option(self, group):
        choices = list(self.widgets[group])
        current_choice = self.vars[group]
        current_idx = choices.index(current_choice)
        next_choice = choices[(current_idx + 1) % len(choices)]

        self.widgets[group][current_choice].choice.set(next_choice)
        self.vars[group] = next_choice


@register_tab
class TabRunAnalysis(Tab):

    name = "run_analysis"

    def init_tab(self):
        cfg = config.there.run_analysis
        if self.vars["rundate"] is None:
            self.vars["rundate"] = date.today() - timedelta(days=1)

        # Banner at top of page
        if cfg.banner.bool:
            self._banner = tk.PhotoImage(file=files.path("there_banner"))
            banner_line = ttk.Frame(self)
            tk.Label(banner_line, image=self._banner).pack(side=tk.TOP, fill=tk.BOTH)
            banner_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Choose necessary options for Where command
        cmd_line = ttk.Frame(self)
        self.add_dropdown(cmd_line, "pipeline_all")
        self.add_dropdown(cmd_line, "rundate")
        self.add_dropdown(cmd_line, "session_all")
        self.add_dropdown(cmd_line, "id_all")
        self.add_checkboxes(cmd_line, None, "force_rerun")
        cmd_line.pack(side=tk.TOP, fill=tk.X, expand=False)

        # Add options
        config_line = ttk.Frame(self)
        cfg_tree = ConfigTree(config_line)
        cfg_tree.update(config.where)
        cfg_tree.pack()
        config_line.pack(side=tk.TOP, fill=tk.X, expand=True)

        # Buttons for running Where
        button_line = ttk.Frame(self)
        self.add_buttons(button_line, self, *cfg.buttons.list)
        button_line.pack(side=tk.BOTTOM, fill=tk.X, expand=False)

        self.update_all()


@register_tab
class TabTimeseries(TabFigure):

    name = "timeseries"

    def init_tab(self):

        # Fix date to 1970-01-01
        rundate = date(1970, 1, 1)
        self.vars.update(dict(rundate=rundate, session="", **config.date_vars(rundate)))

        # Initialize tab as usual
        super().init_tab()


@register_tab
class TabSession(TabFigure):

    name = "session"

    def update_config(self, config_opts):
        cfg_log = [f"Updating configuration:"] + config_opts
        log.info("\n  ".join(cfg_log))

        config_opts = [f"--{o}" for o in config_opts if not o.startswith("--")]
        with config.update_tech_config(use_options=False, **self.vars) as cfg:
            cfg.update_from_options(options=config_opts, source=util.get_program_name())

    def show_log(self):
        log_path = files.path("log", file_vars=self.vars)
        if not log_path.exists():
            log.warn("Log file {} does not exist", log_path)
            return

        log.print_file(log_path, config.there.print_log_level.str.upper())  # TODO: Use log_level enum

    def show_map(self):
        map_path = files.path("output_web_map", file_vars=self.vars)
        if not map_path.exists():
            from where.writers import web_map

            web_map.web_map_writer(self.figure.dataset)
        webbrowser.open(map_path.as_uri())

    def edit_config(self):
        cfg_path = files.path("config", file_vars=self.vars)
        if not cfg_path.exists():
            log.warn(f"Config file '{cfg_path}' does not exist")
            return

        editor.edit(cfg_path)
        setup.add_timestamp(self.vars["rundate"], self.vars["pipeline"], self.vars["session"], "last update")

    button_log = show_log
    button_show_map = show_map
    button_config = edit_config


class UpdateMixin:
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.next_update_func = None

    def set_next_update(self, next_update_func):
        self.next_update_func = next_update_func

    def update_next(self):
        if self.vars["do_update"] and self.next_update_func is not None:
            self.next_update_func()

    @contextmanager
    def no_update(self):
        """Delay updates until end of context manager"""
        self.vars["do_update"] = False
        yield
        self.vars["do_update"] = True


class ConfigTree(ttk.Treeview, UpdateMixin):

    def __init__(self, master):
        super().__init__(master, columns=("value", "help"))
        self.heading("value", text="Value")
        self.column("value", width=300)

        self.heading("help", text="Description")
        self.column("help", width=700)

    def update(self, cfg):
        print("Updating ConfigTree")
        for idx, section in enumerate(cfg.sections):
            tree_section = self.insert("", idx, text=section.name)
            for idx_section, (key, entry) in enumerate(section.items()):
                self.insert(tree_section, idx_section, text=key, values=(entry.str, entry.help))


class Plot(FigureCanvasTkAgg, UpdateMixin):
    def __init__(self, master, vars_):
        figsize = (8, config.there.minimum_figure_height.float)
        self.figure = matplotlib.figure.Figure(figsize=figsize, dpi=100)
        self.graph = self.figure.add_subplot(1, 1, 1)
        super().__init__(self.figure, master)

        self.master = master
        self.filters = list()
        self.vars = vars_
        self.dataset = None
        self.other_dataset = None
        self.next_update_func = None
        self.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        NavigationToolbar2Tk(self, master).update()
        self._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.picker_id = None
        self.ignore_obs = None
        self.cmap = config.there.colormap.str
        self.event_colors = dict()
        self.draw()

    def add_filter(self, filter_):
        self.filters.append(filter_.name.replace("filter_", ""))

    def update_dataset(self):
        if self.vars["rundate"] is None:
            location = config.files.directory_work.directory.replaced.str.split("{")[0]
            log.error(f"No data found in {location}. Run '{where.__executable__}' to generate data.")
            return

        # Read dataset from disk
        self.dataset = data.Dataset(use_options=False, **self.vars)

        # Add event interval field
        events = self.dataset.get_events()
        if events and "time" in self.dataset.fields:
            obs_epochs = self.dataset.time.utc.jd
            event_epochs = np.array([e[0].utc.jd for e in self.dataset.get_events()])
            event, obs = np.meshgrid(event_epochs, obs_epochs)
            if "event_interval" not in self.dataset.fields:
                self.dataset.add_float("event_interval")
            self.dataset.event_interval[:] = np.sum(event - obs < 0, axis=1) + 0.5

        # Observations that will not be plotted
        self.ignore_obs = np.zeros(self.dataset.num_obs, dtype=bool)

        # Update the next widget
        self.update_next()

    def remember_data(self):
        self.other_dataset = self.dataset
        self.vars["x_axis_other"] = self.vars.get("x_axis_data")
        self.vars["y_axis_other"] = self.vars.get("y_axis_data")
        self.vars["color_other"] = self.vars.get("color_data")
        self.vars["size_other"] = self.vars.get("size_data")
        self.update_plot()

    def forget_data(self):
        self.other_dataset = None
        self.vars["x_axis_other"] = None
        self.vars["y_axis_other"] = None
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
        pipeline = self.vars["pipeline"].upper()
        session = f" {self.vars['session']}" if self.vars.get("session") else ""
        date = f" {self.vars['date']}" if self.vars.get("date") else ""
        stage = f" - {self.vars['stage']}".title() if self.vars.get("stage") else ""
        return f"{pipeline}{session}{date}{stage}"

    @property
    def xlabel(self):
        unit = self.dataset.unit(self.vars["x_axis_name"])
        unit_str = " [{}]".format(unit) if unit else ""
        return self.vars["x_axis_name"].replace("_", " ").replace(".", " - ").title() + unit_str

    @property
    def xlim(self):
        return self._calculate_range(self.vars["x_axis_data"])

    @property
    def ylabel(self):
        unit = self.dataset.unit(self.vars["y_axis_name"])
        unit_str = " [{}]".format(unit) if unit else ""
        return self.vars["y_axis_name"].replace("_", " ").replace(".", " - ").title() + unit_str

    @property
    def ylim(self):
        return self._calculate_range(self.vars["y_axis_data"])

    def update_plot(self):
        log.debug("Updating the {}-plot", self.vars["plot_type"])
        tooltip_fields = config.there.tooltip_fields.tuple
        idx, idx_other = self.do_filter()
        self.figure.canvas.mpl_connect("pick_event", self.dbl_click_pick)  # TODO: Delete these?

        # Use the registered plotting functions to plot the correct plot
        self.figure.clear()
        PLOT_TYPES[self.vars["plot_type"]](self)

        # Some quick info, ignore if info can't be calculated
        info_data = self.vars["y_axis_data"]
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
            if event.mouseevent.dblclick:
                return
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
                        self.vars["x_axis_name"],
                        self.vars["y_axis_name"],
                        self.vars["color_name"],
                        self.vars["size_name"],
                    ]
                    if f in self.dataset.fields
                }
            )

            for ind in self.event2dset(event.ind):
                texts = [
                    f"{f}: {self.dataset.plot_values(f)[ind]:{config.FMT_datetime}}"
                    for f in fields
                    if f.startswith("time")
                ]
                texts += ["{}: {}".format(f, self.dataset[f][ind]) for f in fields if not f.startswith("time")]
                log.info("\n       ".join(texts))
                self.master.status("  ".join(texts))

        self.picker_id = self.figure.canvas.mpl_connect("pick_event", on_pick)
        self.figure.subplots_adjust(right=0.99, top=0.95)
        self.draw()

    def do_filter(self):
        filter_dict = {f: self.vars[f"filter_{f}"] for f in self.filters if self.vars[f"filter_{f}"] != "no filter"}
        idx_data = self.dataset.filter(**filter_dict, idx=np.logical_not(self.ignore_obs))
        try:
            idx_other = self.other_dataset.filter(**filter_dict)
        except AttributeError:
            idx_other = idx_data
        return idx_data, idx_other

    def event2dset(self, ind):
        """Convert event.ind indexing to dataset indices

        The event.ind is not counting observations that are filtered out.
        """
        idx, _ = self.do_filter()
        data_vars = ("x_axis_data", "y_axis_data", "color_data", "size_data")
        idx_finite = np.ones(sum(idx), dtype=bool)
        for var in data_vars:
            try:
                data = self.vars[var][idx]
                idx_finite &= np.all(np.isfinite(data), axis=tuple(range(1, data.ndim)))
            except TypeError:
                idx_finite &= np.ones(sum(idx), dtype=bool)

        i_obs = np.where(idx)[0][idx_finite]
        return i_obs[ind]

    def dbl_click_pick(self, event, mouse_event=None):
        if mouse_event is None:
            mouse_event = event.mouseevent

        if not mouse_event.dblclick:
            return False, dict(ind=list())

        dbl_click_func = getattr(self, f"dbl_click_{self.vars['double_click']}")
        return dbl_click_func(event, mouse_event)

    def dbl_click_do_nothing(self, _, __):
        log.debug("Doing nothing about double click")
        return False, dict()

    def dbl_click_add_clock_break(self, event, mouse_event):
        if "ind" in dir(event):  # Workaround so that clock breaks are not added twice
            return False, dict()

        if self.vars["filter_station"] == "no filter":
            log.error("Choose a station to add a clock break")
        else:
            time = Time(mouse_event.xdata, format="plot_date", scale="utc")

            # Check if there is a suspected clock break nearby
            all_data = self.vars["x_axis_data"]
            threshold = (max(all_data) - min(all_data)).total_seconds() * unit.seconds2days / 200
            for suspect_time, _, suspect_station in self.dataset.get_events("suspected_clock_break"):
                if (
                    suspect_station == self.vars["filter_station"]
                    and np.abs(time.utc.jd - suspect_time.utc.jd) < threshold
                ):
                    log.info("Converting suspected clock break")
                    time = suspect_time
                    # TODO: dset.remove_event ...

            clock_break = f"{self.vars['filter_station']} {time.datetime:{config.FMT_datetime}}"
            log.info(f"Adding clock break: '{clock_break}'")

            # Add event on dataset to visualize
            self.dataset.add_event(time, "unprocessed_clock_break", self.vars["filter_station"])
            self.dataset.write()
            self.update_plot()

            # Add to config file
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.vlbi_clock_correction.clock_breaks.as_list(", *")
                updated = ", ".join(sorted(current + [clock_break]))
                cfg.update("vlbi_clock_correction", "clock_breaks", updated, source=util.get_program_name())

        return False, dict()

    def dbl_click_ignore_observation(self, event, _):
        if "ind" not in dir(event):
            return False, dict()

        station_field = "baseline" if "baseline" in self.dataset.fields else "station"

        # Find which observations to ignore
        ignore_list = list()
        for idx in self.event2dset(event.ind):
            ignore_str = f"{self.dataset.time.utc.iso[idx]} {self.dataset[station_field][idx]}"
            log.info(f"Adding '{ignore_str}' to ignore_observation")
            ignore_list.append(ignore_str)

            # Add to ignore filter for visualization
            # TODO: Should ignore_obs be in dset.meta so it will be permanent?
            self.ignore_obs[idx] = True

        self.update_plot()

        # Add to config file
        with config.update_tech_config(use_options=False, **self.vars) as cfg:
            current = cfg.ignore_observation.observations.as_list(", *")
            updated = ", ".join(sorted(current + ignore_list))
            cfg.update("ignore_observation", "observations", updated, source=util.get_program_name())

        return False, dict()

    def dbl_click_go_to_session(self, event, mouse_event):
        if "ind" not in dir(event) or len(event.ind) == 0:
            return False, dict()

        # Figure out info about the observation that was clicked
        idx = self.event2dset(event.ind)[0]  # Use first event index
        var_names = ["rundate", "session", "pipeline", "id"]
        var_names += [f[7:] for f in self.vars.keys() if f.startswith("filter_")]

        session_vars = dict()
        for var in var_names:
            fvar = var[7:] if var.startswith("filter_") else var
            if fvar in self.dataset.fields:
                session_vars[var] = self.dataset[fvar][idx]
            elif fvar in self.vars:
                session_vars[var] = self.vars[fvar]

        # Fix date variables
        if "rundate" in session_vars:
            session_vars["rundate"] = datetime.strptime(session_vars["rundate"], config.FMT_date)
            session_vars["date"] = session_vars["rundate"].strftime("%Y%m%d")
            session_vars["date_session"] = f"{session_vars['date']}/{session_vars['session']}"
        if "id" in session_vars:
            session_vars["id"] = session_vars["id"].lstrip("_")

        # Update session variables
        log.info(f"Opening {session_vars['pipeline'].upper()} {session_vars.get('date_session', '')}")
        self.master.status(f"Opening {session_vars['pipeline'].upper()} {session_vars.get('date_session', '')}")

        session_tab = self.master.master.master.tabs["session"]
        with session_tab.figure.no_update():
            for key, value in session_vars.items():
                if key in session_tab.widgets:
                    session_tab.widgets[key].choice.set(value)
                elif f"filter_{key}" in session_tab.widgets:
                    session_tab.widgets[f"filter_{key}"].choice.set(value)

        # Switch to session tab and update
        session_tab.update_all()
        session_tab.select()

        return False, dict()

    def button_ignore_baseline(self):
        if "filter_baseline" not in self.vars or self.vars["filter_baseline"] == "no filter":
            log.error("Choose a baseline in the filter menu to ignore it")
        else:
            log.info(f"Adding {self.vars['filter_baseline']} to ignore_baseline")
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.vlbi_ignore_baseline.baselines.as_list(", *")
                updated = ", ".join(sorted(current + [self.vars["filter_baseline"]]))
                cfg.update("vlbi_ignore_baseline", "baselines", updated, source=util.get_program_name())

    def button_ignore_station(self):
        if "filter_station" not in self.vars or self.vars["filter_station"] == "no filter":
            log.error("Choose a station in the filter menu to ignore it")
        else:
            log.info(f"Adding {self.vars['filter_station']} to ignore_station")
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.ignore_station.stations.as_list(", *")
                updated = ", ".join(sorted(current + [self.vars["filter_station"]]))
                cfg.update("ignore_station", "stations", updated, source=util.get_program_name())

    def button_ignore_source(self):
        if "filter_source" not in self.vars or self.vars["filter_source"] == "no filter":
            log.error("Choose a source in the filter menu to ignore it")
        else:
            log.info(f"Adding {self.vars['filter_source']} to ignore_source")
            with config.update_tech_config(use_options=False, **self.vars) as cfg:
                current = cfg.vlbi_ignore_source.sources.as_list(", *")
                updated = ", ".join(sorted(current + [self.vars["filter_source"]]))
                cfg.update("vlbi_ignore_source", "sources", updated, source=util.get_program_name())

    # Different types of plots
    @plot_type
    def plot_scatter(self, gs=gridspec.GridSpec(1, 1)[0]):
        idx, idx_other = self.do_filter()
        x_data, x_events = self._identify_events(self.vars["x_axis_data"])
        y_data, y_events = self._identify_events(self.vars["y_axis_data"])

        if x_data.ndim < 1 or y_data.ndim < 1:
            return idx, idx_other

        def event_pick(event, axis):
            all_data = self.vars[axis + "_axis_data"]
            try:
                threshold = (max(all_data) - min(all_data)).total_seconds() * unit.seconds2days / 100
            except AttributeError:
                threshold = (max(all_data) - min(all_data)) / 100
            e_time, e_type, e_description = event
            e_type = e_type.replace("_", " ").title()

            def on_pick(_, mouse_event):
                if mouse_event.dblclick:
                    return False, dict(ind=list())
                mouse_time = getattr(mouse_event, axis + "data")
                if mouse_time and abs(e_time.plot_date - mouse_time) < threshold:
                    log.info(f"Event: {e_time.datetime:{config.FMT_datetime}} - {e_type}\n       {e_description}")
                    self.master.status(f"Event: {e_time.datetime:{config.FMT_datetime}} {e_type}: {e_description}")
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
            xlim = self._pad_range(self.xlim if x_data.ndim == 1 else self.xlim[num_x])
            ylim = self._pad_range(self.ylim if y_data.ndim == 1 else self.ylim[num_y])

            try:
                ax.scatter(
                    self.vars.get("x_axis_other")[idx_other][idx_x],
                    self.vars.get("y_axis_other")[idx_other][idx_y],
                    c=config.there.scatter.color_remember.str,
                    s=self.vars.get("size_other")[idx_other],
                    marker=config.there.scatter.marker_remember.str,
                    alpha=config.there.scatter.alpha_remember.float,
                )
            except (IndexError, TypeError, ValueError):
                log.debug("Not plotting other data")
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
                    (x[0].datetime, x[0].datetime), ylim, ":", color=self.event_colors[x[1]], picker=event_pick(x, "x")
                )
            for y in y_events:
                ax.plot(
                    xlim, (y[0].datetime, y[0].datetime), ":", color=self.event_colors[y[1]], picker=event_pick(y, "y")
                )

            # Label subplot
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
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
        hist_data, hist_range = self._convert_datetime_to_sec(self.vars["x_axis_data"])
        hist_range = self._pad_range(hist_range if hist_data.ndim == 1 else self._flatten_range(hist_range))

        histogram.hist(hist_data[idx], bins=99, range=hist_range, alpha=0.75)
        histogram.set_xlim(hist_range)
        histogram.tick_params(bottom=False, labelbottom=False)
        histogram.set_ylabel("Count")

        # Vertical histogram
        histogram = self.figure.add_subplot(gs[3])
        hist_data, hist_range = self._convert_datetime_to_sec(self.vars["y_axis_data"])
        hist_range = self._pad_range(hist_range if hist_data.ndim == 1 else self._flatten_range(hist_range))

        histogram.hist(hist_data[idx], bins=49, range=hist_range, alpha=0.75, orientation="horizontal")
        histogram.set_ylim(hist_range)
        histogram.tick_params(left=False, labelleft=False)
        histogram.set_xlabel("Count")

    @plot_type
    def plot_polar(self):
        ax = self.figure.add_axes([0.1, 0.05, 0.8, 0.8], polar=True)
        ax.clear()

        idx, idx_other = self.do_filter()
        x_data, x_lim = self._convert_datetime_to_sec(self._project_to_1d(self.vars["x_axis_data"]))
        y_data, y_lim = self._convert_datetime_to_sec(self._project_to_1d(self.vars["y_axis_data"]))

        # Convert radians to degrees for certain fields
        y_name = self.vars["y_axis_name"]
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

        # Pick out data: Prefer y_axis, but alternatively plot x_axis if it is 3d
        if self.vars["y_axis_3d"]:
            plot_data = self.vars["y_axis_data"]
            x_name, y_name, z_name = ["{}: {}".format(xyz, self.ylabel) for xyz in "XYZ"]
        elif self.vars["x_axis_3d"]:
            plot_data = self.vars["x_axis_data"]
            x_name, y_name, z_name = ["{}: {}".format(xyz, self.xlabel) for xyz in "XYZ"]
        else:
            data_x, _ = self._convert_datetime_to_sec(self.vars["x_axis_data"])
            data_y, _ = self._convert_datetime_to_sec(self.vars["y_axis_data"])
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
        plot_range = self._calculate_range(plot_data)

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

    def _calculate_range(self, values):
        """Calculate range of given data
        """
        if self.vars.get("scale_to_filter", False):
            idx_filter, _ = self.do_filter()
        else:
            idx_filter = np.ones(len(values), dtype=bool)

        values = values[idx_filter]
        try:
            idx_finite = np.all(np.isfinite(values), axis=tuple(range(1, values.ndim)))
        except TypeError:  # Datetimes throw a TypeError on isnan/isfinite
            idx_finite = np.ones(len(values), dtype=bool)

        if not np.any(idx_finite):
            return (0, 1) if values.ndim == 1 else ((0, 1),) * values.shape[1]

        if values.ndim == 1:
            return np.min(values[idx_finite]), np.max(values[idx_finite])
        else:
            return tuple(zip(np.min(values[idx_finite], axis=0), np.max(values[idx_finite], axis=0)))

    @staticmethod
    def _flatten_range(ranges):
        """Flatten multidimensional range
        """
        r_min, r_max = zip(*ranges)
        return min(r_min), max(r_max)

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
    def __init__(self, master, name, vars_, figure, **kwargs):
        self.choice = tk.IntVar()
        label = name.replace("_", " ").title()
        super().__init__(master, text=label, command=self.update_vars, variable=self.choice, **kwargs)
        self.name = name
        self.vars = vars_
        self.figure = figure
        self.vars[self.name] = False

    def update_vars(self):
        self.vars[self.name] = bool(self.choice.get())
        if self.figure is not None:
            self.figure.update_plot()
        log.debug("Setting {} = {}", self.name, self.vars[self.name])


class Radiobutton(ttk.Radiobutton):

    groupvars = dict()

    def __init__(self, master, name, group, vars_, **kwargs):
        self.name = name
        self.group = group
        text = name.replace("_", " ").capitalize()
        group = f"{master}.{group}"
        if group not in self.groupvars:
            self.groupvars[group] = tk.StringVar()
            self.groupvars[group].set(self.name)
        self.choice = self.groupvars[group]
        super().__init__(master, text=text, command=self.update_vars, variable=self.choice, value=self.name, **kwargs)
        self.vars = vars_
        self.update_vars()

    def update_vars(self):
        self.vars[self.group] = self.choice.get()
        log.debug("Setting {} = {}", self.group, self.vars[self.group])


class DD_Date(tk.Frame, UpdateMixin):

    name = "no_name"
    width = 0
    start_year = 1971

    def __init__(self, master, vars_, figure, init_val=None, **kwargs):
        super().__init__(master, **kwargs)
        DateVar = namedtuple("DateVar", ["year", "month", "day"])
        self.choices = DateVar(tk.StringVar(name="year"), tk.StringVar(name="month"), tk.StringVar(name="day"))
        self.dropdowns = (
            ttk.Combobox(self, width=5, textvariable=self.choices.year, state="readonly", **kwargs),
            ttk.Combobox(self, width=3, textvariable=self.choices.month, state="readonly", **kwargs),
            ttk.Combobox(self, width=3, textvariable=self.choices.day, state="readonly", **kwargs),
        )
        self.vars = vars_
        self.figure = figure
        self.next_update_func = None

        init_val = self.vars.get(self.var, init_val)
        if not init_val:
            init_strs = ["", "", ""]
        else:
            init_strs = init_val.strftime("%Y-%-m-%-d").split("-")

        for dropdown, var, init in zip(self.dropdowns, self.choices, init_strs):
            dropdown.pack(side=tk.LEFT, padx=3)
            var.set(init.lstrip("0"))
            var.trace("w", self.choose)

        self.update_options()  # Seemingly needed to get correct number of days on startup ...

    @property
    def label(self):
        return self.name.replace("_", " ").title()

    @property
    def var(self):
        return self.name

    def choose(self, name, *_):
        vars_ = self.parse_var(name)
        self.vars.update(vars_)
        if name in ("year", "month"):
            self.update_options()
        for key, value in vars_.items():
            log.debug("Setting {} = {}", key, value)
        self.update_next()

    def parse_var(self, name):
        prefix = f"{self.var}_"
        value = int(getattr(self.choices, name).get())
        value_date = None
        date_vars = {k[len(prefix) :]: v for k, v in self.vars.items() if k.startswith(prefix)}
        date_vars[name] = value
        if "year" in date_vars and "month" in date_vars and "day" in date_vars:
            try:
                value_date = date(date_vars["year"], date_vars["month"], date_vars["day"])
            except ValueError:
                log.warn("{year}-{month}-{day} is not a valid date".format(**date_vars))

        if value_date is None:
            return {f"{prefix}{name}": value}
        else:
            return {self.var: value_date, f"{prefix}{name}": value}

    def update_options(self):
        log.debug("Updating {}", self.name)
        values = (self.read_year(), self.read_month(), self.read_day())
        for dropdown, value, choice in zip(self.dropdowns, values, self.choices):
            dropdown["values"] = value
            if dropdown["values"] and choice.get() not in dropdown["values"]:
                choice.set(dropdown["values"][-1])
            elif not dropdown["values"]:
                choice.set("")
            else:
                choice.set(choice.get())

    def read_year(self):
        return list(range(date.today().year, self.start_year - 1, -1))

    def read_month(self):
        return list(range(1, 12 + 1))

    def read_day(self):
        num_days = 31

        # Figure out number of days based on year and month
        prefix = f"{self.var}_"
        if f"{prefix}year" in self.vars and f"{prefix}month" in self.vars:
            _, num_days = calendar.monthrange(self.vars[f"{prefix}year"], self.vars[f"{prefix}month"])

        return list(range(1, num_days + 1))


@register_dropdown
class DD_RunDate(DD_Date):

    name = "rundate"


class DropdownPicker(ttk.Combobox, UpdateMixin):

    name = "no_name"
    width = 8

    def __init__(self, master, vars_, figure, init_val="", **kwargs):
        self.choice = tk.StringVar()
        super().__init__(master, width=self.width, textvariable=self.choice, state="readonly", **kwargs)

        self.vars = vars_
        self.figure = figure
        self.next_update_func = None
        self.choice.set(self.vars.get(self.var, init_val))
        self.choice.trace("w", self.choose)

    @property
    def label(self):
        return self.name.replace("_", " ").title()

    @property
    def var(self):
        return self.name

    def choose(self, *_):
        vars_ = self.parse_vars()
        self.vars.update(vars_)
        for name, value in vars_.items():
            log.debug("Setting {} = {}", name, value)
        self.update_next()

    def parse_vars(self):
        return {self.var: self.choice.get()}

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


@register_dropdown
class DD_PipelineAll(DropdownPicker):

    name = "pipeline_all"
    label = "Pipeline"
    var = "pipeline"
    width = 8

    def read_options(self):
        """Read possible pipelines
        """
        return pipelines.names()

    def parse_vars(self):
        """TODO: This can be removed when tech is renamed pipeline"""
        return {"tech": self.choice.get(), "pipeline": self.choice.get()}


@register_dropdown
class DD_PipelineDisk(DropdownPicker):

    name = "pipeline"
    width = 8

    def read_options(self):
        """Read possible pipelines from disk, so we only show pipelines that have been analyzed
        """
        return sorted(files.glob_variable("directory_work", "tech", r"[_a-z]+"))

    def parse_vars(self):
        """TODO: This can be removed when tech is renamed pipeline"""
        return {"tech": self.choice.get(), "pipeline": self.choice.get()}


@register_dropdown
class DD_SessionAll(DropdownPicker):

    name = "session_all"
    label = "Session"
    var = "session"
    width = 8

    def read_options(self):
        if not self.vars.get("pipeline") or not self.vars.get("rundate"):
            return list()
        return pipelines.list_sessions(self.vars["rundate"], self.vars["pipeline"])

    def parse_vars(self):
        return {"session": self.choice.get(), "dataset_name": self.choice.get()}


@register_dropdown
class DD_Session(DropdownPicker):

    name = "date_session"
    label = "Session"
    width = 20

    def read_options(self):
        """Read sessions from filenames

        A session in this case is given by a date and a session name.
        """
        file_vars = dict(user=self.vars["user"], tech=self.vars["tech"])
        paths = files.glob_paths("dataset_hdf5", file_vars=file_vars)
        dirs = [p.parent.name for p in paths]
        sessions = set()
        for dirname in dirs:
            session = dirname.split("_")[0]
            if session.startswith("19700101"):
                continue  # Ignore timeseries datasets
            sessions.add(f"{session[:8]}/{session[8:]}")
        return sorted(sessions)

    def parse_vars(self):
        """Construct date variables and split out session name
        """
        date_session = self.choice.get()
        date_value, _, session = date_session.partition("/")
        rundate = datetime.strptime(date_value, "%Y%m%d").date() if date_value else None
        session_vars = config.date_vars(rundate)
        session_vars["rundate"] = rundate
        session_vars["date"] = date_value
        session_vars["date_session"] = date_session
        session_vars["session"] = session
        return session_vars


@register_dropdown
class DD_IdAll(DropdownPicker):

    name = "id_all"
    label = "Id"
    var = "id"
    width = 25

    def read_options(self):
        """Read ids from config"""
        ids = [""] + config.there.run_analysis.ids.list
        return sorted(i.lstrip("_") for i in ids)

    def parse_vars(self):
        """Add leading underscore to id
        """
        id_ = self.choice.get()
        return {self.var: f"_{id_}" if id_ else ""}


@register_dropdown
class DD_Id(DropdownPicker):

    name = "id"
    width = 25

    def read_options(self):
        """Read ids from filenames

        """
        file_vars = {k: v for k, v in self.vars.items() if k not in ("id", "stage", "dataset_name", "dataset_id")}
        ids = files.glob_variable("dataset_hdf5", "id", r"|_[_\w]+", file_vars=file_vars)
        return sorted(i.lstrip("_") for i in ids)

    def parse_vars(self):
        """Add leading underscore to id
        """
        id_ = self.choice.get()
        return {self.name: f"_{id_}" if id_ else ""}


@register_dropdown
class DD_Stage(DropdownPicker):

    name = "stage"
    width = 14

    def read_options(self):
        """Read pipeline and stage from filenames
        """
        file_vars = {k: v for k, v in self.vars.items() if k not in ("stage",)}
        stages = files.glob_variable("dataset_hdf5", "stage", r"[a-z_]+", file_vars=file_vars)
        return sorted(stages, key=self._sorter)

    def _sorter(self, stage):
        """Return a sort key for the given stage

        Sorts first on tech then on stage such that later stages are sorted first. Unknown stages are sorted
        alphabetically at the end. Stages can be prioritized by adding them to stages in there.conf.

        Args:
            tech_stage:  Tuple with tech and stage.

        Returns:
            Tuple that can be sorted on.
        """
        tech = self.vars["tech"]
        stages = config.there.session.stages.list + pipelines.stages(tech)[::-1]  # Reversed to sort later stages first
        try:
            stage_id = stages.index(stage)
        except ValueError:
            stage_id = len(stages)  # Unknown stages are sorted last
        return (stage_id, stage)


@register_dropdown
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


@register_dropdown
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

        return {self.name + "_data": values, self.name + "_name": field, self.name + "_3d": do_3d_plot}


@register_dropdown
class DD_XAxis(DD_Field):

    name = "x_axis"
    label = "X-Axis"


@register_dropdown
class DD_YAxis(DD_Field):

    name = "y_axis"
    label = "Y-Axis"


@register_dropdown
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


@register_dropdown
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
        name = f"filter_{filter_name}"
        label = filter_name.replace("_", " ").title()

        def __init__(self, master, file_vars, figure, **kwargs):
            super().__init__(master, file_vars, figure, init_val="no filter", **kwargs)
            figure.add_filter(self)

        def read_options(self):
            try:
                return ["no filter"] + self.figure.dataset.unique(self.name[7:])
            except KeyError:
                return ["no filter"]

    register_dropdown(DD_Filter)
    return DD_Filter


def normalize(values, figure):
    """Normalize the values of values so they are in the [0, 1]-interval
    """
    if figure.vars.get("scale_to_filter", False):
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
