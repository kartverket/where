import os

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

from vlbi_latlon import latlon

SIMPLE_MAP = True
ANIMATION = False
DPI = 100

# Read names and coordinates
lats = [v[0] for v in latlon.values()]
lons = [v[1] for v in latlon.values()]
names = [k for k in latlon.keys()]

idx_1 = names.index("Ny-Alesund")
idx_2 = names.index("Hartebeesthoek")

# import IPython; IPython.embed()
plt.xkcd()


if SIMPLE_MAP:
    # Add a simple map of all VLBI stations with names
    map = Basemap(projection="robin", lon_0=0, resolution="c")  # robin, moll
    map.bluemarble()
    map.drawmapboundary()
    #    map.drawmeridians(np.arange(0, 360, 30))
    #    map.drawparallels(np.arange(-90, 90, 30))

    x, y = map(lons, lats)
    for xp, yp, name in zip(x, y, names):
        print(xp, yp, name)
        if 16000000 < xp < 21000000 and 10000000 < yp < 15000000:
            continue
        plt.text(
            xp - 500000,
            yp,
            name,
            color="red",
            verticalalignment="center",
            horizontalalignment="right",
            fontsize=6,
            weight="bold",
        )
    map.plot(x, y, "ro")

    plt.savefig(os.path.join("figure", "vlbi_kart.png"), bbox_inches="tight", transparent=True, dpi=1.5 * DPI)
    plt.close()


if ANIMATION:
    for lon_0 in range(0, 360, 10):
        filename = os.path.join("figure", "vlbi_baseline_{:03d}{}.png".format(lon_0, "{}"))
        print(filename)
        map = Basemap(projection="ortho", lat_0=20, lon_0=lon_0, resolution="l", area_thresh=1000)
        map.bluemarble()

        # draw the edge of the map projection region (the projection limb)
        map.drawmapboundary()
        # draw lat/lon grid lines every 30 degrees.
        map.drawmeridians(np.arange(0, 360, 30))
        map.drawparallels(np.arange(-90, 90, 30))

        # compute the native map projection coordinates for stations.
        x, y = map(lons, lats)
        map.plot(x, y, "ro", alpha=0.5)

        # plot lines towards a radio source for the two given stations
        x1, x2 = x[idx_1], x[idx_2]
        y1, y2 = y[idx_1], y[idx_2]

        map.plot(x1, y1, "ro")
        map.plot(x2, y2, "ro")
        plt.savefig(filename.format("_pts"), bbox_inches="tight", transparent=True, dpi=DPI)

        # At least one of the two given stations are on the back side of the earth
        if np.max((x1, x2, y1, y2)) > 1e15:
            plt.close()
            continue

        map.plot((0, x1), (0, y1), "r-")
        map.plot((x2 - x1, x2), (y2 - y1, y2), "r-")
        plt.savefig(filename.format("_src"), bbox_inches="tight", transparent=True, dpi=DPI)

        map.drawgreatcircle(lons[idx_1], lats[idx_1], lons[idx_2], lats[idx_2], color="r", alpha=0.5)
        map.plot((x1, x2), (y1, y2), "r--", linewidth=2)
        plt.savefig(filename.format("_baseline"), bbox_inches="tight", transparent=True, dpi=DPI)

        # normal vector
        len_station_1 = np.sqrt(x1 ** 2 + y1 ** 2)
        len_baseline = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        angle = (
            np.arccos((x1 * (x1 - x2) + y1 * (y1 - y2)) / (len_station_1 * len_baseline)) * np.sign(x2 - y2 * x1 / y1)
        )
        len_normal = len_baseline * np.sin(angle)
        xn = x2 - y1 * len_normal / len_station_1
        yn = y2 + x1 * len_normal / len_station_1

        map.plot((x1, xn), (y1, yn), "r-", linewidth=2)
        map.plot((x2, xn), (y2, yn), "r-", linewidth=2)

        #        import IPython; IPython.embed()
        plt.savefig(filename.format(""), bbox_inches="tight", transparent=True, dpi=DPI)
        plt.close()
