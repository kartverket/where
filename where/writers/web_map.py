"""Create a map showing all sites in analysis and some simple statistics

Description:
------------

asdf

"""

# Standard library imports

# External library imports
from branca import colormap
import folium
import numpy as np

# Where imports
from where.lib import files
from where.lib import log
from where.lib import plugins
from where.lib.unit import unit


@plugins.register
def web_map_writer(dset):
    file_path = files.path("output_web_map", file_vars=dset.vars)
    log.info("Storing a web map at '{}'. Open in a browser to look at it", file_path)

    sites = read_site_latlons(dset)
    map = draw_map(dset, sites)
    map.save(str(file_path))


def read_site_latlons(dset):
    sites = dict()
    for station in dset.unique("station"):
        for _ in dset.for_each("station"):
            llh = dset.first("site_pos.llh", station=station)
            if llh is not None:
                sites[station] = (llh[0] * unit.rad2degree, llh[1] * unit.rad2degree)

    return sites


def draw_map(dset, latlons):
    stations = dset.unique("station")
    rms = dset.rms("residual") * unit.m2mm

    # Map layers
    map = folium.Map((20, 0), zoom_start=2)  # Default is OpenStreetMap
    folium.TileLayer("Stamen Toner").add_to(map)
    folium.TileLayer("CartoDB Positron").add_to(map)

    # Colors
    colors = colormap.LinearColormap(("green", "yellow", "red"), vmin=rms / 2, vmax=rms * 2)
    colors.caption = "RMS [mm]"
    map.add_child(colors)

    # Stations
    stations_layer = folium.FeatureGroup(name="Stations")
    for sta in stations:
        idx = dset.filter(station=sta)
        rms = dset.rms("residual", idx=idx) * unit.m2mm
        other_idx = [dset.filter(idx=idx, station=other) for other in stations]
        others = [
            "{:.3f} mm to {} ({})".format(dset.rms("residual", idx=i) * unit.m2mm, other, sum(i))
            for i, other in zip(other_idx, stations)
            if other != sta and any(i)
        ]
        popup = "<b>{}</b> ({}): {:.3f} mm<br />{}".format(sta, sum(idx), rms, "<br />".join(others))
        stations_layer.add_child(
            folium.CircleMarker(
                latlons[sta],
                popup=popup,
                fill=True,
                color=colors(rms),
                radius=max(5, 10 * np.sqrt(np.mean(idx) * len(stations))),
            )
        )
    map.add_child(stations_layer)

    # Baselines
    for sta_1 in stations:
        baseline_layer = folium.FeatureGroup(name="{} baselines".format(sta_1))
        for sta_2 in stations:
            idx = dset.filter(station=sta_1) & dset.filter(station=sta_2)
            if sta_1 == sta_2 or not any(idx):
                continue
            rms = dset.rms("residual", idx=idx) * unit.m2mm
            locations = [latlons[s] for s in (sta_1, sta_2)]
            popup = "<b>{} - {}</b> ({}): {:.3f} mm".format(sta_1, sta_2, sum(idx), rms)
            baseline_layer.add_child(
                folium.PolyLine(
                    locations, popup=popup, color=colors(rms), weight=max(1, 2 * np.mean(idx) * len(stations) ** 2)
                )
            )
        if baseline_layer._children:
            map.add_child(baseline_layer)

    folium.LayerControl().add_to(map)
    return map
