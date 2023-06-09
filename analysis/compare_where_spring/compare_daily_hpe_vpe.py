# Standard library imports
from datetime import datetime
import pathlib

# External library imports
import numpy as np
import pandas as pd

# Midgard imports
from midgard.plot.matplotlib_extension import plot

FIGURE_FORMAT = "png"


def _get_diff(x_arrays1, y_arrays1, x_arrays2, y_arrays2, stations):
    """Determine difference
    """
    y_arrays = list()
    for station, y1, y2 in zip(stations, y_arrays1, y_arrays2):
        diff = y1 - y2
        y_arrays.append(diff)
        print(
            f"{station.upper()}: mean = {np.nanmean(diff):5.2f}, min = {np.nanmin(diff):5.2f}, max = {np.nanmax(diff):5.2f}, rms = {np.sqrt(np.nanmean(np.square(diff))):5.2f}"
        )

    return x_arrays1, y_arrays


def _get_spring_pickle_data(paths, solution, frequency, stations):
    """Get x- and y-arrays for plotting SPRING solution
    """

    # Merge different monthly solutions together
    df_merged = pd.DataFrame()
    for path in paths:
        if df_merged.empty:
            df_merged = pd.read_pickle(path)
        else:
            df_merged = pd.concat([df_merged, pd.read_pickle(path)])

    if frequency == "e1_nequick":
        df = df_merged[["Epoch", "Station", solution]][df_merged["Type"] == "SPV_NEQ_L1_GAL"]
    elif frequency == "e1e5a":
        df = df_merged[["Epoch", "Station", solution]][df_merged["Type"] == "rih_GAL_DF_E1_5A_POS"]

    x_arrays = list()
    y_arrays = list()
    for station in stations:
        idx = df["Station"] == station.upper()
        date = np.array([datetime.strptime(epoch, "%d/%m/%Y") for epoch in df[idx]["Epoch"]])
        x_arrays.append(date)
        print(df[idx])  # DEBUG
        y_arrays.append(np.array(df[idx][solution]))

    return x_arrays, y_arrays


def _get_table_data(path, stations):
    """Get x- and y-arrays for plotting based on data table
    """
    df = pd.read_table(path, delim_whitespace=True)

    x_arrays = list()
    y_arrays = list()
    date = np.array([datetime.strptime(epoch, "%d-%m-%Y") for epoch in df["date"]])
    for station in stations:
        x_arrays.append(date)
        data = np.array(df[station.lower()])
        y_arrays.append(data)
        print(
            f"{station.upper()}: mean = {np.nanmean(data):5.2f}, min = {np.nanmin(data):5.2f}, max = {np.nanmax(data):5.2f}, rms = {np.sqrt(np.nanmean(np.square(data))):5.2f}"
        )

    return x_arrays, y_arrays


def _plot(x_arrays, y_arrays, stations, field, solution, figure_path) -> None:
    """Generate plots
    """
    opt_args = {
        "colormap": "tab20",
        "figsize": (7, 3),
        # "grid": True,
        "marker": "o",
        "markersize": "3",
        "linestyle": "solid",
        "plot_to": "file",
        "plot_type": "plot",
        # "statistic": ["rms", "mean", "std", "min", "max", "percentile"], #TODO: Is only shown for data, which are plotted at last.
        # "title": title,
        #"xlim": [datetime(2019, 7, 1), datetime(2019, 9, 30)],
    }

    if field == "hpe":
        opt_args.update({"ylim": [0.0, 8.0]})
    elif field == "vpe":
        opt_args.update({"ylim": [0.0, 10.0]})
    elif field == "vel_3d":
        opt_args.update({"ylim": [0.0, 0.09]})

    if solution == "Difference":
        if field in ["hpe", "vpe"]:
            opt_args.update({"ylim": [-0.5, 1.0]})
        elif field == "vel_3d":
            opt_args.update({"ylim": [-0.04, 0.0]})

    #colors = ["gold", "orange", "red", "violet", "blue", "green"]
    colors = ["orange", "red", "violet", "blue", "green"]

    # Generate plot
    plot(
        x_arrays=x_arrays,
        y_arrays=y_arrays,
        xlabel="Time [GPS]",
        ylabel=f"{solution}: 3D VE 95%" if field == "vel_3d" else f"{solution}: {field.upper()} 95%",
        y_unit="m",
        labels=stations,
        colors=colors,
        figure_path=figure_path,
        opt_args=opt_args,
    )


############################################################################################
#
#
#                                MAIN PROGRAM
#
#
############################################################################################
if __name__ == "__main__":

    stations = ["nabd", "hons", "vegs", "hofs", "krss"]
    directory = pathlib.Path("/home/dahmic/where/analysis/compare_where_spring")
    paths = [
        # directory / 'data/spring_pickle/2019_07_spring_e1_nequick_position.dat',
        # directory / 'data/spring_pickle/2019_08_spring_e1_nequick_position.dat',
        # directory / 'data/spring_pickle/2019_09_spring_e1_nequick_position.dat',
        #directory / "data/spring_pickle/2019_07_spring_e1e5a_position.dat",
        directory / 'data/spring_pickle/2020_01_spring_e1_nequick_velocity.dat',
    ]
    compare_against = "spring"  # spring, terrapos
    sol = "e1_nequick"  # e1_nequick, e1e5a, e1e5b

    field_def = {
#        "hpe": {
#            "where": directory / f"data/2019_07_09_where_{sol}_hpe.dat",
#            #'where': directory / f"data/2019_07_where_{sol}_hpe.dat",
#            "spring": directory / f"data/2019_07_09_spring_{sol}_hpe.dat",
#            "spring_pickle": "2D 95pct",
#            #'terrapos': directory / f"data/2019_07_09_terrapos_e1_klob_hpe.dat",
#            "terrapos": directory / f"data/2019_07_09_terrapos_{sol}_hpe.dat",
#        },
#        "vpe": {
#            "where": directory / f"data/2019_07_09_where_{sol}_vpe.dat",
#            #'where': directory / f"data/2019_07_where_{sol}_hpe.dat",
#            "spring": directory / f"data/2019_07_09_spring_{sol}_vpe.dat",
#            "spring_pickle": "Vertical 95pct",
#            #'terrapos': directory / f"data/2019_07_09_terrapos_e1_klob_vpe.dat",
#            "terrapos": directory / f"data/2019_07_09_terrapos_{sol}_vpe.dat",
#        },
        "vel_3d": {
            "where": directory / f"data/2020_01_03_where_{sol}_vel_3d.dat",
            "spring": directory / f"data/2020_01_03_spring_{sol}_vel_3d.dat",
            #"spring_pickle": "Vertical 95pct",
        },
    }

    for field in field_def.keys():

        print(f"Plot {field.upper()} solution.")

        # Read SPRING data
        # x_s, y_s = _get_spring_pickle_data(paths, field_def[field]['spring_pickle'], sol, stations)

        # Read Where data
        print(f"Read Where solution: {field_def[field]['where']}")
        x_w, y_w = _get_table_data(field_def[field]["where"], stations)

        if compare_against == "spring":
            # Read Spring data
            print(f"Read SPRING solution: {field_def[field]['spring']}")
            x_comp, y_comp = _get_table_data(field_def[field]["spring"], stations)

        if compare_against == "terrapos":
            # Read Terrapos data
            print(f"Read TERRAPOS solution: {field_def[field]['terrapos']}")
            x_comp, y_comp = _get_table_data(field_def[field]["terrapos"], stations)

        # Determine difference
        print("Generate differenced solution.")
        x_diff, y_diff = _get_diff(x_w, y_w, x_comp, y_comp, stations)

        # Plot results
        if compare_against == "spring":
            _plot(
                x_comp,
                y_comp,
                stations,
                field,
                "Spring",
                directory / f"plots/plot_{field}_spring_{sol}.{FIGURE_FORMAT}",
            )
        if compare_against == "terrapos":
            _plot(
                x_comp,
                y_comp,
                stations,
                field,
                "Terrapos",
                directory / f"plots/plot_{field}_terrapos_{sol}.{FIGURE_FORMAT}",
            )
            # _plot(x_comp, y_comp, stations, field, 'Terrapos', directory/f"plots/plot_{field}_terrapos_e1_klob.{FIGURE_FORMAT}")

        _plot(x_w, y_w, stations, field, "Where", directory / f"plots/plot_{field}_where_{sol}.{FIGURE_FORMAT}")
        _plot(
            x_diff,
            y_diff,
            stations,
            field,
            "Difference",
            directory / f"plots/plot_{field}_where_{compare_against}_{sol}.{FIGURE_FORMAT}",
        )
