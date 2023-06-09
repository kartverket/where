#!/usr/bin/env python3
"""Compare Where and gLAB analysis concerning difference between precise and broadcast orbits

Description:
------------

NOTE: gLAB uses simplified expressions for determination of solar coordinates (based on Section 3.3.2 in
      :cite:`montenbruck2012`), whereas Where uses more sophisticated JPL ephemeris. Also the broadcast velocities are
      determined based on the difference between satellite positions before and after the needed epoch. In addition
      gLAB corrects for Earth rotation, which leads to large velocity differences and based on that to differences in
      the yaw reference coordinate system (e.g. orb_diff_acr.itrs -> difference on meter level for along-track and
      cross-track between Where and gLAB solution). Why is gLAB correcting for Earth rotation?

1. Run gLAB SISRE analysis
cd ~/SOURCE/GNSS/gLAB_3.0.0_Linux
./gLAB_linux -input:nav ~/where/data/obs/orb/brdc/2016/brdm0610.16p
           -input:SP3 ~/where/data/obs/orb/igs/1886/com18862.sp3
           -input:ant ~/where/data/apriori/antenna/igs14.atx
           -pre:dec 300 -model:clock:deg 1 -model:orbit:deg 10 
           > ~/where/src/analysis/compare_where_glab/20160301_orbdiff_glab.orig

./gLAB_linux -input:nav ~/where/data/obs/orb/brdc/2016/stas0610.16l
           -input:SP3 ~/where/data/obs/orb/igs/1886/com18862.sp3
           -input:ant ~/where/data/apriori/antenna/igs14.atx
           -pre:dec 300 -model:clock:deg 1 -model:orbit:deg 10 
           > ~/where/src/analysis/compare_where_glab/20160301_orbdiff_glab_stas_E.orig

#Navigation file filtered by Where
./gLAB_linux -input:nav ~/where/data/obs/orb/brdc/2016/BRDW00IGS_R_20160610000_01D_MN.rnx
           -input:SP3 ~/where/data/obs/orb/igs/1886/com18862.sp3
           -input:ant ~/where/data/apriori/antenna/igs14.atx
           -pre:dec 300 -model:clock:deg 1 -model:orbit:deg 10 
           > ~/where/src/analysis/compare_where_glab/20160301_orbdiff_glab.orig

cd ~/where/src/analysis/compare_where_glab
grep SATDIFF 20160301_orbdiff_glab.orig > 20160301_orbdiff_glab_satdiff.out
grep 'BRDC:' 20160301_orbdiff_glab.orig > 20160301_orbdiff_glab_brdc.out

2. Run Where SISRE analysis

where 2016 3 1 --sisre

GENERAL CONFIGURATION:
sampling_rate                  = 300
systems                        = G E
satellites                     = G01 G02 G03 G04 G05 G06 G07 G08 G09 G10 G11
                                 G12 G13 G14 G15 G16 G17 G18 G19 G20 G21 G22
                                 G23 G24 G25 G26 G27 G28 G29 G30 G31 G32
use_mixed_brdc_file            = True
brdc_block_nearest_to          = transmission_time_positive

CONFIGURATION files.conf:
[gnss_rinex_nav_M]
filename        = brdm{$doy}0.{$yy}p{$gz}
directory       = {$path_data}/obs/orb/brdc/{$yyyy}

[gnss_orbit_sp3]
filename        = com{$gpsweek}{$dow}.sp3
directory       = {$path_data}/obs/orb/igs/{$gpsweek}

[gnss_antex]
filename        = igs14.atx
directory       = {$path_data}/apriori/antenna


3. Read gLAB output file and Where Dataset

4. Generate gLAB Dataset

5. Compare Where and gLAB Dataset

Authors:
--------
* Michael Daehnn <michael.daehnn@kartverket.no>

"""
from sys import exit
import collections
from datetime import date, datetime, timedelta
from astropy.time import TimeDelta
import numpy as np
from where import data


def _difference_where_glab(dwhere, dglab):

    # Get common dataset
    hash_glab = [(t.gps.isot + sn) for t, sn in dglab.values("time", "satellite")]
    hash_where = [(t.gps.isot + sn) for t, sn in dwhere.values("time", "satellite")]
    hash_both = sorted(set(hash_glab) & set(hash_where))
    idx_w = [hash_where.index(h) for h in hash_both]
    idx_g = [hash_glab.index(h) for h in hash_both]

    ddiff = data.Dataset(date, "sisre", "orbdiff_where_gLAB", "stas", 0, empty=True)
    ddiff.num_obs = len(hash_both)

    common_fields = ["sisre", "sisre_orb", "orb_diff_3d", "clk_diff"]

    for field in common_fields:
        if field in dwhere.fields:
            if field in dglab.fields:
                diff = dwhere[field][idx_w] - dglab[field][idx_g]
                ddiff.add_float(field, val=diff)

    ddiff.add_time("time", val=dwhere.time.gps[idx_w], scale="gps")
    ddiff.add_position("orb_diff", time="time", itrs=dwhere.orb_diff.itrs[idx_w] - dglab.orb_diff.itrs[idx_g])
    ddiff.add_position(
        "orb_diff_acr", time="time", itrs=dwhere.orb_diff_acr.itrs[idx_w] - dglab.orb_diff_acr.itrs[idx_g]
    )
    ddiff.add_position(
        "diff_sat_phase_center",
        time="time",
        itrs=dwhere.gnss_satellite_phase_center_offset.itrs[idx_w]
        - dglab.gnss_satellite_phase_center_offset.itrs[idx_g],
    )
    ddiff.add_float("satnum", val=dglab.satnum[idx_g])
    ddiff.add_text("satellite", val=np.array(dglab.satellite)[idx_g])
    ddiff.add_text("system", val=np.array(dglab.system)[idx_g])
    ddiff.write()


def _generate_glab_dset(inpath, date):
    """Generate gLAB dataset based on SATDIFF output (function compareOrbits() in gLAB.c )

    Read gLAB output file.
    """
    sys_ids = list()
    sys_id_definition = {"BDS": "C", "GAL": "E", "GPS": "G"}

    # Read gLAB output
    float_fields = dict(
        year=1,
        doy=2,
        seconds=3,
        satnum=5,
        sisre=6,
        sisre_orb=7,
        orb_diff_3d=8,
        clk_diff=9,
        orb_diff_radial=10,
        orb_diff_along_track=11,
        orb_diff_cross_track=12,
        orb_diff_x=13,
        orb_diff_y=14,
        orb_diff_z=15,
        sat_phase_center_offset_x=16,
        sat_phase_center_offset_y=17,
        sat_phase_center_offset_z=18,
    )
    system_field = dict(system=4)
    fields = float_fields.copy()
    fields.update(system_field)

    # glab_data = np.loadtxt(inpath, usecols=fields.values())
    glab_data = np.genfromtxt(inpath, usecols=fields.values(), dtype="str")

    # Define gLAB dataset
    dglab = data.Dataset(date, "sisre", "gLAB", "stas", 0, empty=True)
    dglab.num_obs = glab_data.shape[0]

    for idx, f in enumerate(float_fields.keys()):
        dglab.add_float(f, val=glab_data[:, idx])

    for system in glab_data[:, 17]:
        try:
            sys_ids.append(sys_id_definition[system])
        except KeyError:
            print("Fatal error: System identifier '{}' is not defined.".format(system))
            exit()
    dglab.add_text("system", val=sys_ids)

    # Calculated fields
    epochs = list()
    for idx, (year, doy, seconds) in enumerate(dglab.values("year", "doy", "seconds")):
        epochs.append(datetime.strptime("{:.0f} {:.0f}".format(year, doy), "%Y %j") + timedelta(seconds=seconds))

    dglab.add_time("time", val=epochs, format="datetime", scale="gps")
    dglab.add_position(
        "orb_diff_acr",
        time="time",
        itrs=np.vstack((dglab.orb_diff_along_track, dglab.orb_diff_cross_track, dglab.orb_diff_radial)).T,
    )
    dglab.add_position(
        "orb_diff", time="time", itrs=np.vstack((dglab.orb_diff_x, dglab.orb_diff_y, dglab.orb_diff_z)).T
    )
    dglab.add_position(
        "gnss_satellite_phase_center_offset",
        time="time",
        itrs=np.vstack(
            (dglab.sat_phase_center_offset_x, dglab.sat_phase_center_offset_y, dglab.sat_phase_center_offset_z)
        ).T,
    )

    # Get satellite identifier
    satnum = dglab.satnum.astype(int)
    satnum = np.core.defchararray.zfill(satnum.astype(str), 2)
    satellites = np.core.defchararray.add(dglab.system, satnum)
    dglab.add_text("satellite", val=satellites)
    dglab.write()

    return dglab


############################################################################################
#
#
#                                MAIN PROGRAM
#
#
############################################################################################
if __name__ == "__main__":
    # year = 2010
    # month = 1
    # day = 25
    year = 2016
    month = 3
    day = 1
    date = date(year, month, day)

    # Read gLAB output file
    inpath = "{:%Y%m%d}_orbdiff_glab_satdiff.out".format(date)
    # inpath = '{:%Y%m%d}_orbdiff_glab_satdiff_stas_E.out'.format(date)

    # Where dataset
    dwhere = data.Dataset(date, "sisre", "orb_diff", "stas", 0)

    # gLAB dataset
    dglab = _generate_glab_dset(inpath, date)

    # Determine difference between the Where and gLAB dataset
    _difference_where_glab(dwhere, dglab)
