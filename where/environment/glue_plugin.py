"""Package for handling Where data in Glue

Functions for plugging Where data into the Glue vizualization tool. To use the plugin, add the following line to your
~/.glue/config.py-file:

    from where.environment import glue_plugin

If you do not have glue installed, see http://glueviz.org/ for instructions. If you have glue, but not the
~/.glue/config.py file, run the program `glue-config`.


"""

from datetime import datetime
import os.path
import re
import uuid

from qtpy import QtWidgets
from glue.config import data_factory, menubar_plugin
from glue.core import Data
from glue.utils.qt import update_combobox

from where import data
from where.lib import config
from where.lib import files
from where.lib import time


def is_where_dataset(filename, **kwargs):
    """Check if the given file is a Where dataset

    For now we do a simple check on the filename ending. This should be made more robust.

    @todo More robust check for whether a given file belongs to a Where dataset.

    Args:
        filename:  Given filename.

    Returns:
        Bool, True if filename is part of a Where dataset, False otherwise.
    """
    return filename.endswith(".json") or filename.endswith("hdf5")


@data_factory("Where Dataset", is_where_dataset, priority=10000)
def read_where_dataset_from_file(filename):
    """Read a Where dataset from a filename and return Glue Data-objects

    This function is used when Glue opens a file that is a Where dataset.

    Args:
        filename:   Filename (either .json or .hdf5) representing a Where dataset.

    Returns:
        List of Glue Data-objects, one for each Where dataset in the given file.
    """
    print("-> read_where_dataset(", filename, ")")

    dsets = data.list_datasets_from_filename(filename)
    return [open_where_dataset_as_glue(d) for d in dsets]


def open_where_dataset_as_glue(dataset_vars):

    dset = data.Dataset(**dataset_vars)

    # Set where config
    import where

    rundate = dataset_vars["rundate"]
    where.set_config(rundate.year, rundate.month, rundate.day, dataset_vars["tech"])

    # Add fields of dataset as components to glue
    components = dict()

    # As 1-d tables
    for field in dset.fields:
        try:
            values = dset[field]
        except:
            print("Not able to add", field)
            continue

        if isinstance(values, time.Time):
            values = values.mjd

        try:
            suffixes = "x y z" if values.shape[1] == 3 else "1_x 1_y 1_z 2_x 2_y 2_z"
            for column, suffix in zip(values.T, suffixes.split()):
                components[field + "_" + suffix] = column
                print("Adding {}".format(field + "_" + suffix))
        except (AttributeError, IndexError):
            try:
                len(values)
            except TypeError:
                print("Not able to add", field)
                continue

            components[field] = values
            print("Adding {}".format(field))

    glue_data = Data(**components)
    glue_data.label = "{tech} {stage} {rundate:%Y%m%d} {name}".format(name=dset.name, **dataset_vars)
    return glue_data


@menubar_plugin("Where toolbar")
def load_where_dataset(session, data_collection):
    print("-> load_where_dataset(", session, data_collection, ")")
    selector = DataSelectorWidget(data_collection)
    toolbar = QtWidgets.QToolBar()
    toolbar.addWidget(selector)
    session.application.addToolBar(toolbar)


def simple_update_combobox(combo, labels):
    """
    Wrapper around glue's update_combobox, which automatically sets the
    userData to a unique ID.
    """
    print("-> simple_update_combobox(", combo, labels, ")")
    labels_with_data = [(label, str(uuid.uuid4())) for label in labels]
    return update_combobox(combo, labels_with_data)


class UpdateMixin():

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.next_update_func = None

    def set_next_update(self, next_update_func):
        self.next_update_func = next_update_func

    def update_next(self):
        if self.next_update_func is not None:
            self.next_update_func()


class DataSelectorWidget(QtWidgets.QWidget):

    def __init__(self, data_collection=None, parent=None):
        print("-> DataSelectorWidget.__init__(", data_collection, parent, ")")
        super().__init__(parent=parent)

        # Keep reference to data collection
        self.data_collection = data_collection
        self.dset_vars = dict(location="work")
        self.widget_updates = list()

        # Set up the layout and widgets
        self.setup_ui()

        # We trigger the changed method to populate the combo boxes
        self.update_all()

    def add_dropdown(self, layout, dropdown):
        print("-> DataSelectorWidget.add_dropdown(", layout, dropdown, ")")
        layout.addWidget(QtWidgets.QLabel(dropdown.name.title() + ":"))
        layout.addWidget(dropdown)
        layout.addStretch()
        self.add_widget_update(dropdown, dropdown.update)

    def add_widget_update(self, widget, update_func):
        if self.widget_updates:
            last_widget = self.widget_updates[-1][0]
            last_widget.set_next_update(update_func)
        self.widget_updates.append((widget, update_func))

    def update_all(self):
        self.widget_updates[0][1]()

    def setup_ui(self):
        print("-> DataSelectorWidget.setup_ui()")

        layout = QtWidgets.QHBoxLayout()

        # Add dropdowns
        self.add_dropdown(layout, DD_Location(self, self.dset_vars))
        self.add_dropdown(layout, DD_User(self, self.dset_vars))
        self.add_dropdown(layout, DD_Date(self, self.dset_vars))
        self.add_dropdown(layout, DD_Technique(self, self.dset_vars))
        self.add_dropdown(layout, DD_Dataset(self, self.dset_vars))

        # Add button
        update_button = QtWidgets.QPushButton("Update")
        layout.addWidget(update_button)
        update_button.clicked.connect(self.update_where_dataset_from_toolbar)

        layout.setContentsMargins(5, 5, 5, 5)
        self.setLayout(layout)

    def update_where_dataset_from_toolbar(self):
        glue_data = open_where_dataset_as_glue(self.dset_vars)
        glue_data._is_where_plugin_data = True

        glue_data_to_update = None
        for data in self.data_collection:
            if hasattr(data, "_is_where_plugin_data") and data._is_where_plugin_data:
                glue_data_to_update = data
                break

        if glue_data_to_update is None:
            self.data_collection.append(glue_data)
        else:
            glue_data_to_update.update_values_from_data(glue_data)


class DropdownPicker(QtWidgets.QComboBox, UpdateMixin):

    name = "no_name"
    width = 8

    def __init__(self, parent, dset_vars):
        print("-> DropdownPicker(", parent, dset_vars, ")")
        super().__init__(parent=parent)
        self.dset_vars = dset_vars
        self.model_vars = config.files.vars
        self.currentIndexChanged.connect(self.choose)

    def choose(self):
        vars_ = self.changed()
        self.dset_vars.update(vars_)
        for name, value in vars_.items():
            print("Setting {} = {}".format(name, value))
        self.update_next()

    def changed(self):
        print("-> DropdownPicker.changed()", self.name)
        return {self.name: self.currentText()}

    def update(self):
        simple_update_combobox(self, list())


class DD_Location(DropdownPicker):

    name = "location"
    width = 8

    def update(self):
        simple_update_combobox(self, sorted([d[10:] for d in config.files.sections if d.startswith("directory_")]))


class DD_User(DropdownPicker):

    name = "user"
    width = 8

    def update(self):
        """Read users from the file directories
        """
        users = files.glob_variable("dataset_hdf5", "user", r"[a-z]+")
        simple_update_combobox(self, sorted(users))


class DD_Date(DropdownPicker):

    name = "date"
    width = 30

    def update(self):
        """Read dates from filenames
        """
        vars_ = dict(self.model_vars, user=self.dset_vars["user"])
        paths = files.glob_paths("dataset_hdf5", file_vars=vars_)
        dirs = [os.path.basename(os.path.dirname(p)) for p in paths]
        dates = set()
        for dirname in dirs:
            dates.add(dirname.replace("_", "/", 1))

        simple_update_combobox(self, sorted(dates))

    def changed(self):
        """Construct date variables and split out timestamp
        """
        parts = (self.currentText() + "/").split("/")
        rundate = datetime.strptime(parts[0], "%Y%m%d") if parts[0] else None
        vars_ = config.date_vars(rundate) if rundate else dict()
        vars_["rundate"] = rundate
        vars_["date"] = parts[0]
        vars_["timestamp"] = parts[1].split("_")[0]
        try:
            # Test if timestamp is valid
            datetime.strptime(vars_["timestamp"], config.FMT_dt_file)
            parts[1] = "_".join(parts[1].split("_")[1:])
        except ValueError:
            vars_["timestamp"] = ""
        vars_["id"] = "_" + parts[1] if parts[1] else ""
        return vars_


class DD_Technique(DropdownPicker):

    name = "technique"
    width = 14

    def update(self):
        """Read technique and stage from filenames
        """
        vars_ = {k: v for k, v in self.dset_vars.items() if k not in ("tech", "stage")}
        paths = files.glob_paths("dataset_hdf5", file_vars=vars_)
        simple_update_combobox(
            self, ["{0}/{2}".format(*s) for s in [os.path.basename(p).split("-") for p in sorted(paths)]]
        )

    def changed(self):
        if not self.currentText():
            return dict()

        tech, stage = self.currentText().split("/")
        return dict(tech=tech, stage=stage)


class DD_Dataset(DropdownPicker):

    name = "dataset"
    width = 14
    next_update_func = None

    def update(self):
        """Read dataset name and id
        """
        print("-> DD_Dataset.update()", self.dset_vars)
        simple_update_combobox(self, sorted(data.list_datasets(**self.dset_vars), reverse=True))

    def changed(self):
        dataset_name, dataset_id = self.currentText().split("/")
        return dict(dataset_name=dataset_name, dataset_id=int(dataset_id))
