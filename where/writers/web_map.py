"""Create a web map showing all sites in analysis and some simple statistics

Description:
------------

The web map is based on folium and shows the position of each site in the dataset.
The colors in the plot indicate the RMS of the individual stations and baselines.
The line width in the plot indicate the number of observations on the baseline.

"""

# Standard library imports

# External library imports
from branca import colormap
import folium
import numpy as np

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit

# Where imports
from where.lib import config
from where.lib import log


@plugins.register
def web_map_writer(dset):
    file_vars = {**dset.vars, **dset.analysis}
    file_path = config.files.path("output_web_map", file_vars=file_vars)
    log.info(f"Storing a web map at '{file_path}'. Open in a browser to look at it")

    sites = read_site_latlons(dset)
    webmap = draw_map(dset, sites)
    webmap.save(str(file_path))


def read_site_latlons(dset):
    sites = dict()
    for station in dset.unique("station"):
        sites[station] = (dset.meta[station]["latitude"] * Unit.rad2deg, dset.meta[station]["longitude"] * Unit.rad2deg)
    return sites


def draw_map(dset, latlons):
    stations = dset.unique("station")
    rms = dset.rms("residual") * Unit.m2mm

    # Map layers
    webmap = folium.Map((20, 0), zoom_start=2)  # Default is OpenStreetMap
    folium.TileLayer("Stamen Toner").add_to(webmap)
    folium.TileLayer("CartoDB Positron").add_to(webmap)

    # Colors
    colors = colormap.LinearColormap(("green", "yellow", "red"), vmin=rms / 2, vmax=rms * 2)
    colors.caption = "RMS [mm]"
    webmap.add_child(colors)

    # Stations
    stations_layer = folium.FeatureGroup(name="Stations")
    for sta in stations:
        idx = dset.filter(station=sta)
        rms = dset.rms("residual", idx=idx) * Unit.m2mm
        other_idx = [dset.filter(idx=idx, station=other) for other in stations]
        others = [
            f"{dset.rms('residual', idx=i) * Unit.m2mm:.3f} mm to {other} ({sum(i)})"
            for i, other in zip(other_idx, stations)
            if other != sta and any(i)
        ]
        popup = f"<b>{sta}</b> ({sum(idx)}): {rms:.3f} mm<br />{'<br />'.join(others)}"
        stations_layer.add_child(
            folium.CircleMarker(
                latlons[sta],
                popup=popup,
                fill=True,
                color=colors(rms),
                radius=max(5, 10 * np.sqrt(np.mean(idx) * len(stations))),
            )
        )
    webmap.add_child(stations_layer)

    # Baselines
    for sta_1 in stations:
        baseline_layer = folium.FeatureGroup(name="{sta_1} baselines")
        for sta_2 in stations:
            idx = dset.filter(station=sta_1) & dset.filter(station=sta_2)
            if sta_1 == sta_2 or not any(idx):
                continue
            rms = dset.rms("residual", idx=idx) * Unit.m2mm
            locations = [latlons[s] for s in (sta_1, sta_2)]
            popup = "<b>{sta_1} - {sta_2}</b> ({sum(idx)}): {rms:.3f} mm"
            baseline_layer.add_child(
                folium.PolyLine(
                    locations, popup=popup, color=colors(rms), weight=max(1, 2 * np.mean(idx) * len(stations) ** 2)
                )
            )
        if baseline_layer._children:
            webmap.add_child(baseline_layer)

    folium.LayerControl().add_to(webmap)
    return webmap
