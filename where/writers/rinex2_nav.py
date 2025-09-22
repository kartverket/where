"""Write RINEX navigation file

Description:
------------
Write data in the Rinex navigation file format 2.11 (see :cite:`rinex2`).

TODO: DOES NOT WORK


"""
from datetime import timedelta

# Midgard imports
from midgard.dev import plugins

# Where imports
from where import apriori
from where.lib import config, log


@plugins.register
def rinex2_nav(dset):
    """Write RINEX navigation file for 3 days (rundate +/-1 day)

    Args:
        dset:       Dataset, a dataset containing the data.
    """

    log.fatal("'rinex2_nav' writer does not work at the moment. Update your configuration and start process without "
              "using this writer.")
    

    # Overwrite Dataset. This is necessary if the writer is called from a analysis (e.g. SISRE) with does not include
    # broadcast ephemeris information.
    # TODO: Is that the best solution?
    if dset.vars["pipeline"] != "rinex_nav":
        brdc = apriori.get(
            "orbit",
            rundate=dset.analysis["rundate"],
            system=tuple(dset.unique("system")),
            station=dset.vars["station"],
            apriori_orbit="broadcast",
        )
        meta = brdc.dset_edit.meta[dset.analysis["rundate"].strftime("%Y-%m-%d")]
        data = brdc.dset_edit  # TODO: Another possibility: brdc.dset_raw

    else:
        meta = dset.meta[dset.analysis["rundate"].strftime("%Y-%m-%d")]
        data = dset  # TODO: Another possibility: brdc.dset_raw

    file_path = config.files.path("output_rinex2_nav", file_vars={**dset.vars, **dset.analysis})
    log.info(f"Write file {file_path}.")

    with config.files.open_path(file_path, create_dirs=True, mode="wt") as fid:

        #
        # Write RINEX navigation header
        #
        if meta["file_type"] == "N":
            file_type = "NAVIGATION DATA"

        fid.write("{:>9s}{:11s}{:40s}RINEX VERSION / TYPE\n".format(meta["version"], "", file_type))
        fid.write(
            "{:20s}{:20s}{:20s}PGM / RUN BY / DATE\n".format(meta["program"], meta["run_by"], meta["file_created"])
        )

        if "comment" in meta:
            for line in meta["comment"]:
                fid.write("{:60s}COMMENT\n".format(line))
        fid.write(
            "{:>14.4e}{:>12.4e}{:>12.4e}{:>12.4e}{:10s}ION ALPHA\n"
            "".format(
                meta["iono_para"]["GPSA"]["para"][0],
                meta["iono_para"]["GPSA"]["para"][1],
                meta["iono_para"]["GPSA"]["para"][2],
                meta["iono_para"]["GPSA"]["para"][3],
                "",
            )
        )
        fid.write(
            "{:>14.4e}{:>12.4e}{:>12.4e}{:>12.4e}{:10s}ION BETA\n"
            "".format(
                meta["iono_para"]["GPSB"]["para"][0],
                meta["iono_para"]["GPSB"]["para"][1],
                meta["iono_para"]["GPSB"]["para"][2],
                meta["iono_para"]["GPSB"]["para"][3],
                "",
            )
        )
        # TODO fid.write('{:>22.12e}{:>19.12e}{:>9d}{:>9d}{:1s}DELTA-UTC: A0,A1,T,W\n'
        #          ''.format(meta['a0'], meta['a1'], int(meta['t']), int(meta['w']), ''))
        fid.write("{:>6d}{:54s}LEAP SECONDS\n".format(int(meta["leap_seconds"]["leap_seconds"]), ""))
        fid.write("{:60s}END OF HEADER\n".format(""))

        #
        # Write RINEX navigation data
        #
        # TODO:
        #        for drow in data.get_rows():
        #            fid.write('{d.sat:2d} {d.time:%Y%m%d %H%M%S} {d.inc0:13.4f}'.format(d=drow))
        #            fid.write('  {d.hurra:14.10f} ...'.format(d=drow))

        for idx in range(0, data.num_obs):
            sat = int(data.satellite[idx][1:3])
            d = data.time.gps.datetime[idx]
            fid.write(
                "{:2d}{:>3s}{:>3d}{:>3d}{:>3d}{:>3d}{:>5.1f}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(
                    sat,
                    str(d.year)[2:4],
                    d.month,
                    d.day,
                    d.hour,
                    d.minute,
                    d.second,
                    data.sat_clock_bias[idx],
                    data.sat_clock_drift[idx],
                    data.sat_clock_drift_rate[idx],
                )
            )
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.iode[idx], data.crs[idx], data.delta_n[idx], data.m0[idx])
            )
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.cuc[idx], data.e[idx], data.cus[idx], data.sqrt_a[idx])
            )
            # TODO: toe depends on GNSS system time -> for BeiDou it has to be changed
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.toe.gps.gpssec[idx], data.cic[idx], data.Omega[idx], data.cis[idx])
            )
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.i0[idx], data.crc[idx], data.omega[idx], data.Omega_dot[idx])
            )
            # TODO: gnss_week depends on GNSS -> for BeiDou it has to be changed
            # TODO: codes_l2 only valid for GPS and QZSS -> Galileo data_source; rest None
            # TODO: 'G': 'l2p_flag', 'J': 'l2p_flag'
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.idot[idx], data.codes_l2[idx], data.gnss_week[idx], data.l2p_flag[idx])
            )
            # TODO: 'G': 'iodc', 'J': 'iodc', 'E': 'bgd_e1_e5b', 'C': 'tgd_b2_b3'
            # TODO: 'G': 'tgd', 'J': 'tgd', 'E': 'bgd_e1_e5a', 'C': 'tgd_b1_b3', 'I': 'tgd'
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.sv_accuracy[idx], data.sv_health[idx], data.tgd[idx], data.iodc[idx])
            )
            # TODO: transmission_time depends on GNSS system time -> for BeiDou it has to be changed
            # TODO: fit_interval only valid for GPS and QZSS -> for BeiDou age_of_clock_corr; rest None
            fid.write(
                "{:22.12e}{:>19.12e}{:>19.12e}{:>19.12e}\n"
                "".format(data.transmission_time.gps.gpssec[idx], data.fit_interval[idx], 0.0, 0.0)
            )
