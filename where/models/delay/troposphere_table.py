"""Calculate total zenith delay

Description:
------------

This model determines the total zenith delay from an external source (e.g. GNSS, ray-tracing).

References:
-----------

[1] Mendes, V.B. and E.C. Pavlis, 2004,
    "High-accuracy zenith delay prediction at optical wavelengths,"
    Geophysical Res. Lett., 31, L14602, doi:10.1029/2004GL020308, 2004

[2] Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
    IERS Technical Note No. 36, BKG (2010)



"""
# External library imports
import numpy as np
from midgard.math import interpolation

# Where imports
from where.ext import iers_2010 as iers
from where.lib import plugins
from where.lib.unit import unit
from where.models.delay.troposphere_radio import gmf_mapping_function, gpt2w, apg_gradient_model, saastamoinen_zenith_hydrostatic_delay, saastamoinen_zenith_wet_delay
from where import parsers
from where.lib import log
from where.lib.time import Time

@plugins.register
def troposphere_for_all_stations(dset):
    """Calculate tropospheric delay for all stations

    Depending on the used technique a Dataset can include one or more stations. In the case of VLBI a baseline between
    two stations are analysed. Here the relative tropospheric delay is needed between two stations and not the absolute
    one.  Therefore if a Dataset has two stations then the difference between the tropospheric delay is calculated by
    using ``multiplier``, which is '-1' for station 1 and '1' for station 2. For one or more than two stations in a
    Dataset the ``multiplier`` is equal to '1'.

    Args:
        dset (Dataset):     Model data.

    Returns:
        numpy.ndarray: Corrections in meters for each observation.
    """
    dset_out = np.zeros(dset.num_obs)
    for multiplier in dset.for_each("station"):
        dset_out += multiplier * troposphere(dset)

    return dset_out

def troposphere(dset):

    # Parse SINEX TRO file
    trop_data = parsers.parse_key("trop_files").as_dict()
    pos = {k: v["pos"] for k, v in trop_data.items()}
    dist = {k: np.linalg.norm(p - dset.site_pos.itrs_pos, axis=1) for k, p in pos.items()}
    trop_site=[min(dist, key=lambda x: dist[x][idx]) for idx in range(dset.num_obs)]
    site_dist=[dist[s][i] for i,s in enumerate(trop_site) ]
    
    pressure, temperature, _, tm, e, _, _, lambd, _ = gpt2w(dset)
    latitude, _, height = dset.site_pos.llh.T
    mh, mw = gmf_mapping_function(dset)
    # Note: Be aware, that the apg.f function uses c = 0.0031 instead of c = 0.0032.
    mg = 1 / (np.sin(dset.site_pos.elevation) * np.tan(dset.site_pos.elevation) + 0.0032)
    zhd = saastamoinen_zenith_hydrostatic_delay(pressure, latitude, height)
    zwd = saastamoinen_zenith_wet_delay(latitude, height, temperature, e)

    gn, ge = apg_gradient_model(dset)
    
    for idx, (site,dist) in enumerate(zip(trop_site,site_dist)):
        if dist > 1000:
            log.warn(f"No GNSS station found within 1000 m of {dset.station[idx]}")
            continue
        sec=Time(trop_data[trop_site[idx]]["epoch"],scale="gps").sec_of_day
        ge[idx] = interpolation.interpolate(sec, trop_data[trop_site[idx]]["tgetot"]/1000, dset.time.gps.sec_of_day, kind="lagrange")[idx]
        gn[idx] = interpolation.interpolate(sec, trop_data[trop_site[idx]]["tgntot"]/1000, dset.time.gps.sec_of_day, kind="lagrange")[idx]
        
    # Line of sight delay
    dT = mh * zhd + mw * zwd + mg * (gn * np.cos(dset.site_pos.azimuth) + ge * np.sin(dset.site_pos.azimuth))
    
    terms_and_levels = dict(
        dT="detail",
        mh="detail",
        zhd="detail",
        mw="operational",
        zwd="detail",
        mg="operational",
        gn="detail",
        ge="detail",
    )
    for term, write_level in terms_and_levels.items():
        field = "troposphere_{}{}".format(term, (dset.default_field_suffix or ""))
        if field in dset.fields:
            dset[field][:] = locals()[term]
        else:
            dset.add_float(field, table="troposphere", val=locals()[term], write_level=write_level)
    #import IPython
    #IPython.embed()
    return dT