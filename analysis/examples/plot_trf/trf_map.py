#!/usr/bin/env python3
"""Plot station coordinates for one reference frame for several techniques

Usage:

    ./trf_map.py tech1 [tech2 ...] [filter1 filter2 ...] [--trf=...]
                 [--no-labels] [--web [--no-clusters]] [-h|--help]

By default, the time series are plotted for the {trf}-reference frame.

One or more of the techniques {techs} should be specified at the command
line. Arguments that are not techniques will be interpreted as filters so that
only stations matching these filters will be plotted on the map.

Additionally, the following options are recognized:

-h, --help      - Show this help text.
--trf=...       - Set which reference frame to use.
--no-labels     - Do not add labels to the plot.
--web           - Save as web map, instead of showing a static map.
--no-clusters   - Do not use clustering in web map.


Installation of GeoPandas and Folium:
--------------------------

We use GeoPandas to do the plotting. Unfortunately, there seems to be some
missing dependencies in the conda version of fiona, which is a dependency of
geopandas. The following seems to work, though:

    conda install -c conda-forge geopandas icu=58 kealib=1.4.7 matplotlib=2

Folium is used to create the web map. It should be straight forward to install
through conda-forge:

    conda install -c conda-forge folium


Examples:
---------

Plot a map of the VLBI network

    ./trf_map.py vlbi

Plot the GNSS network without labels

    ./trf_map.py gnss --no-labels

Plot ITRF for all techniques

    ./trf_map.py vlbi slr gnss doris

Make a web map of the SLR network

    ./trf_map.py slr --web

Make a web map of all techniques without clustering

    ./trf_map.py vlbi slr gnss doris --web --no-clusters

Show stations matching `73` and `ny` on a web map

    ./trf_map.py vlbi slr gnss doris --web 73 ny


Authors:
--------

* Geir Arne Hjelle <geir.arne.hjelle@kartverket.no>

"""

# Standard library imports
from datetime import datetime
import html
import itertools
import sys

# External library imports
from fiona import crs
from shapely import geometry
import geopandas as gpd
import numpy as np
import matplotlib.pyplot as plt

# Where imports
from where import apriori
from where.lib import config
from where.lib import log
from where.ext import sofa
from where import techniques
from where.data.time import Time

MARKERSIZES = iter(range(30, 0, -5))
PLACEMENTS = iter((("left", "center"), ("right", "center"), ("center", "top"), ("center", "bottom")))
COLORS = iter(("red", "blue", "green", "purple"))


def read_trf(tech, trf_names, time, all_pos):
    """Read positions for all sites for the technique and add to dict"""
    markersize = next(MARKERSIZES)
    ha, va = next(PLACEMENTS)
    color = next(COLORS)

    trf = apriori.get("trf", time=time, reference_frames=trf_names)
    for site in trf:
        mean_pos = site.pos.mean(axis=0)
        lon, lat, height = mean_pos.llh
        all_pos.setdefault("tech", list()).append(tech)
        all_pos.setdefault("site_code", list()).append(site.key)
        all_pos.setdefault("name", list()).append(site.name)
        all_pos.setdefault("source", list()).append(site.source)
        all_pos.setdefault("x", list()).append(mean_pos.x)
        all_pos.setdefault("y", list()).append(mean_pos.y)
        all_pos.setdefault("z", list()).append(mean_pos.z)
        all_pos.setdefault("lat", list()).append(np.degrees(lat))
        all_pos.setdefault("lon", list()).append(np.degrees(lon))
        all_pos.setdefault("height", list()).append(height)
        all_pos.setdefault("markersize", list()).append(markersize)
        all_pos.setdefault("ha", list()).append(ha)
        all_pos.setdefault("va", list()).append(va)
        all_pos.setdefault("color", list()).append(color)


def plot_map_static(df_pos, labels=True):
    """Plot a static map using matplotlib"""
    world = gpd.read_file(gpd.datasets.get_path("naturalearth_lowres"))
    ax = world.plot(color="lightgreen", edgecolor="green")
    df_pos.plot(ax=ax, column="tech", legend=True, markersize=df_pos.markersize, linewidth=0)
    if labels:
        df_pos.apply(lambda x: ax.annotate(s=x.site_code, xy=x.geometry.coords[0], ha=x.ha, va=x.va), axis=1)
    plt.show()


def plot_map_web(df_pos, clusters=True):
    """Plot a web map using folium"""
    import folium
    import folium.plugins

    # Add map and map tiles
    map = folium.Map([20, 0], zoom_start=2)  # OpenStreetMap included by default
    folium.TileLayer(
        "http://opencache.statkart.no/gatekeeper/gk/gk.open_gmaps?layers=topo4&zoom={z}&x={x}&y={y}",
        attr='&copy; <a href="http://www.kartverket.no/Kart/Gratis-kartdata/Lisens/">Kartverket</a>',
        name="kartverket",
    ).add_to(map)
    folium.TileLayer("Stamen Terrain").add_to(map)
    folium.TileLayer("Stamen Toner").add_to(map)
    folium.TileLayer("Stamen Watercolor").add_to(map)
    folium.TileLayer("CartoDB Dark_matter").add_to(map)
    folium.TileLayer("CartoDB Positron").add_to(map)

    for tech in df_pos.tech.unique():
        idx = df_pos.tech == tech
        log.info("Adding {} {} sites to the map", sum(idx), tech)

        locations, popups, colors = [], [], []
        for _, row in df_pos[idx].iterrows():
            locations.append((row.geometry.y, row.geometry.x))
            popups.append(
                "<b>{r.site_code} ({r.tech})</b><br>(x, y, z) = ({r.x:.2f}, {r.y:.2f}, {r.z:.2f})<br>"
                "{description}".format(r=row, description=html.escape(row["name"]))
            )
            colors.append(row.color)

        layer = folium.FeatureGroup(name=tech)
        if clusters:
            icons = [folium.map.Icon(color=c, icon="rss", prefix="fa") for c in colors]
            layer.add_child(folium.plugins.MarkerCluster(locations=locations, popups=popups, icons=icons))
        else:
            for location, popup, color in zip(locations, popups, colors):
                layer.add_child(folium.CircleMarker(location, popup=popup, fill=True, color=color))

        map.add_child(layer)
    folium.LayerControl().add_to(map)

    file_path = "trf_map.html"
    log.info("Saving web map to HTML-file {}", file_path)
    map.save(file_path)


def _split_args(args=None):
    """Split command line arguments in positional and optional arguments"""
    if args is None:
        args = sys.argv[1:]

    pos_args, opt_args = list(), list()
    for arg in args:
        if arg.startswith("-"):
            opt_args.append(arg)
        else:
            pos_args.append(arg)

    return pos_args, opt_args


def main():
    log.init()
    trf = "itrf:2014"
    pos_args, opt_args = _split_args()

    # Help
    if len(pos_args) < 1 or "-h" in opt_args or "--help" in opt_args:
        print(__doc__.format(trf=trf, techs=", ".join(techniques.names())))
        sys.exit()

    # Techniques
    techs = sorted(set(pos_args) & set(techniques.names()))

    # TRF
    trfs = [o[6:] for o in opt_args if o.startswith("--trf=")]
    if trfs:
        trf = trfs.pop()
    log.info("Using TRF {}", trf)

    # Use time equal beginning of last year, unless otherwise told
    years = [float(o[7:]) for o in opt_args if o.startswith("--year=")] + [datetime.now().year - 1]
    time = Time([years[0]], fmt="decimalyear", scale="utc")
    log.info("Using time epoch {}", time[0].datetime)

    # Read TRF for each technique
    all_pos = dict()
    for tech in techs:
        log.info("Reading {} for {} from file", trf, tech.upper())
        # Set technique (this should be improved ...)
        config.set_analysis(rundate=None, tech=tech)
        config.files.profiles = [tech]
        read_trf(tech, trf, time, all_pos)

    # Convert to GeoDataFrame
    all_pos["geometry"] = [geometry.Point(lon, lat) for lat, lon in zip(all_pos["lat"], all_pos["lon"])]
    df_pos = gpd.GeoDataFrame(all_pos, geometry="geometry")
    df_pos.crs = crs.from_epsg(4326)

    # Filter
    filters = set(pos_args) - set(techs)
    idx = np.zeros(len(df_pos), dtype=bool)
    if filters:
        for f, field in itertools.product(filters, ("name", "site_code")):
            idx |= df_pos[field].str.lower().str.contains(f.lower())
        log.info("Filtering on '{}' leaves {} sites", ", ".join(filters), sum(idx))
    if not any(idx):
        idx = np.ones(len(df_pos), dtype=bool)

    # Plot map
    log.info("Mapping {} for {}".format(trf, ", ".join(t.upper() for t in techs)))
    if "--web" in opt_args:
        plot_map_web(df_pos[idx], clusters="--no-clusters" not in opt_args)
    else:
        plot_map_static(df_pos[idx], labels="--no-labels" not in opt_args)


if __name__ == "__main__":
    main()
