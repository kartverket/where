
import argparse
import matplotlib
import matplotlib.pyplot as plt
from datetime import date
import itertools
from matplotlib import gridspec
import numpy as np
import sys

import midgard

from where.data import dataset3 as dataset
from where.lib import config

matplotlib.rcParams["axes.formatter.useoffset"] = False


parser = argparse.ArgumentParser(epilog=f"Exapmle: python {sys.argv[0]} 2009 11 2 R1403 time.utc.mjd site_pos.trs --station NYALES20")
parser.add_argument("yyyy", help="Four digit year of session", type=int)
parser.add_argument("mm", help="Month of session", type=int)
parser.add_argument("dd", help="Day of session", type=int)
parser.add_argument("session_code", help="Code of session", type=str)
parser.add_argument("x_field", help="Name of x-axis field (without numbered suffix)", type=str)
parser.add_argument("y_field", help="Name of y-axis field (without numbered suffix)", type=str) 
parser.add_argument("--station", help="Name of station", type=str, default="RAEGSMAR")
parser.add_argument("--id", help="Dataset id of result files.", type=str, default="")
parser.add_argument("--stage", help="Dataset stage", type=str, default="calculate")
parser.add_argument("--label", help="Dataset label ", type=str, default="last")
parser.add_argument("--save_fig", help="Enable this flag to save the plot to file instead of getting an interactive plot", action="store_true", default=False)
args = parser.parse_args()

yyyy = args.yyyy
mm = args.mm
dd = args.dd
session_code = args.session_code
station = args.station
x_field = args.x_field
y_field = args.y_field


stage = args.stage
label = args.label
dset_id = args.id
save_fig = args.save_fig


def insert_suffix(dset, field, suffix):
    field_list = field.split(".")
    f = field_list[0]
    if f not in dset.fields:
        field_list[0] = f"{f}_{suffix}" 
    field = ".".join(field_list)
    return field

def _pad_range(val_range, factor=0.01):
    """Pad value range at both sides
    """
    rmin, rmax = val_range
    try:
        delta = np.fmax((rmax - rmin) * factor, 1e-6)
    except TypeError:  # np.fmax fails on timedeltas
        delta = (rmax - rmin) * factor

    return rmin - delta, rmax + delta


rundate = date(yyyy, mm, dd)

# Setup config
pipeline = "vlbi"
config.set_analysis(rundate=rundate, pipeline=pipeline)
config.files.profiles = [pipeline]
config.read_pipeline(pipeline)

dset = dataset.Dataset.read(rundate=rundate, pipeline=pipeline, stage=stage, label=label, session_code=session_code, id=dset_id)

y_field_1 = insert_suffix(dset, y_field, "1")
y_field_2 = insert_suffix(dset, y_field, "2")
x_field_1 = insert_suffix(dset, x_field, "1")
x_field_2 = insert_suffix(dset, x_field, "2")

idx_1 = dset.filter(station_1=station)
idx_2 = dset.filter(station_2=station)
x_1 = dset[x_field_1][idx_1]
y_1 = dset[y_field_1][idx_1]
x_2 = dset[x_field_2][idx_2]
y_2 = dset[y_field_2][idx_2]
x = np.concatenate((x_1, x_2))
y = np.concatenate((y_1, y_2))

try:
    x_unit = dset.unit_short(x_field_1)
except midgard.dev.exceptions.UnitError:
    x_unit = ""

try:
    y_unit = dset.unit_short(y_field_1)
except midgard.dev.exceptions.UnitError:
    y_unit = ""

footnote = f"source: {config.files.path('dataset', {**dset.vars, **dset.analysis}).name}"     
print(footnote)

if len(x) == 0:
    print(f"No data for {x_field} for {station} in session")
    sys.exit(1)

if len(y) == 0:
    print(f"No data for {y_field} for {station} in session")
    sys.exit(1)


# Handle multi-dimensional data as separate plots
fig = plt.figure()
gs = gridspec.GridSpec(1, 1)[0]
ncols = 1 if x.ndim <= 1 else x.shape[1]
nrows = 1 if y.ndim <= 1 else y.shape[1]
sub_gs = gridspec.GridSpecFromSubplotSpec(nrows, ncols, subplot_spec=gs)

for plot_num, (num_y, num_x) in enumerate(itertools.product(range(nrows), range(ncols))):
    ax = fig.add_subplot(sub_gs[plot_num])
    idx_x = slice(None) if x.ndim == 1 else (slice(None), num_x)
    idx_y = slice(None) if y.ndim == 1 else (slice(None), num_y)
    ax.xaxis.get_major_formatter().set_scientific(False)
    ax.yaxis.get_major_formatter().set_scientific(False)
    ax.scatter(x[idx_x], y[idx_y], marker=".") 
    ax.set_xlim(_pad_range((min(x[idx_x]), max(x[idx_x]))))
    ax.set_ylim(_pad_range((min(y[idx_y]), max(y[idx_y]))))

    if nrows > 1:
        try:
            names = dset[y_field_1].column_names
            ax.set_ylabel(names[num_y])
        except:
            pass

    if ncols > 1:
        try:
            names = dset[x_field_1].column_names
            ax.set_xlabel(names[num_x])
        except:
            pass

fig.text(0.5, 0.01, f"{x_field} {x_unit}", ha="center")
fig.text(0.01, 0.5, f"{y_field} {y_unit}", va="center", rotation="vertical")
fig.align_ylabels()
fig.suptitle(station)    
fig.tight_layout()

if save_fig:
    plt.savefig(f"{station}_{x_field}_{y_field}_{yyyy}{mm}{dd}{session}_{stage}_{label}_{dset_id}.png", dpi=150)
else:
    plt.show()

