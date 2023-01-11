"""Add network volume to dataset

Description:
------------

Computes the network volume and saves the parameter in meta.

Network volume is a useful proxy to say something about how well a session is suited to estimate EOPs. See: 
    "On comparison of the Earth orientation parameters obtained from different VLBI networks and observing programs" 
    by Z. Malkin
    J Geod (2009) 83:547–556
    DOI 10.1007/s00190-008-0265-2
"""
# External library imports
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.spatial import Delaunay

# Midgard imports
from midgard.dev import plugins
from midgard.math.unit import Unit
from midgard.math.constant import constant

# Where imports
from where import apriori
from where.lib import config
from where.lib import log

# Name of section in configuration
_SECTION = "_".join(__name__.split(".")[-1:])


@plugins.register
def vlbi_network_volume(dset: "Dataset") -> None:
    """Add network volume to dataset

    Args:
        dset:     A Dataset containing model data.
    """
    reference_frames = config.tech.reference_frames.str
    trf = apriori.get("trf", time=dset.time.utc.mean, reference_frames=reference_frames)
    site_ids = dset.unique("site_id")
    
    if len(site_ids) < 4:
        log.warn(f"Unable to compute network volume. Only {len(site_ids)} usable stations in session.")
        dset.meta.add("value", np.nan, section="network_volume")
        dset.meta.add("__unit__", "Megameter**3", section="network_volume")
        return
    
    pos = np.zeros((len(site_ids), 3))
    for i, site_id in enumerate(site_ids):    
        pos[i, :] = trf[site_id].pos.trs.val
    
    triang = Delaunay(pos)
    volume = 0

    for t in triang.simplices:
        a, b, c, d = triang.points[t]
        # Volume of tetrahedron #https://en.wikipedia.org/wiki/Tetrahedron
        volume += np.abs(np.dot(a -d, np.cross(b - d, c - d)))/6

    volume = volume * Unit.meter2Megameter**3
    
    dset.meta.add("value", volume, section="network_volume")
    dset.meta.add("__unit__", "Megameter**3", section="network_volume")
    log.info(f"Network volume: {volume:.4f} Mm**3")
    
    if config.tech.vlbi_network_volume.plot_triangulation.bool:
        plot_triangulation(dset, triang, volume)
        
def plot_triangulation(dset, triang, volume):
    """Plot the network polyhedron 
    
    Plots a wireframe earth and the polyhedron resulting from the delaunay triangulation. The plot is saved as a file.
    """
    
    unit_text = r"[$\mathregular{Mm^3}$]"
    title = f"Volume: {volume:12.2f} {unit_text}"
    
    fig = plt.figure(figsize=(12,12), dpi=150)
    ax = fig.gca(projection="3d")
    for t in triang.simplices:
        x, y, z = triang.points[t].T
        ax.plot_trisurf(x,y,z)
        ax.scatter(x,y,z, color="black")

    # Plot earth wireframe as a sphere
    r = constant.a_E # Earth radius 
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x = r * np.cos(u)*np.sin(v)
    y = r * np.sin(u)*np.sin(v)
    z = r * np.cos(v)
    ax.plot_wireframe(x, y, z, color="r", alpha=0.2)
    ax.set_xlabel("meter")
    ax.set_ylabel("meter")
    ax.set_zlabel("meter")
    plt.title(title)

    filename = config.files.path("vlbi_network_triangulation", file_vars=dset.vars)
    if not filename.parent.exists(): 
        os.makedirs(filename.parent)
    plt.savefig(filename, bbox_inches='tight')
    plt.close()




            
